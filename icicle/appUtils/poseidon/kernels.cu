#include "poseidon.cuh"

namespace poseidon {
  template <typename S, int T>
  __global__ void prepare_poseidon_states(S* states, size_t number_of_states, S domain_tag, bool aligned)
  {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int state_number = idx / T;
    if (state_number >= number_of_states) { return; }
    int element_number = idx % T;

    S prepared_element;

    // Domain separation
    if (element_number == 0) {
      prepared_element = domain_tag;
    } else {
      if (aligned) {
        prepared_element = states[idx];
      } else {
        prepared_element = states[idx - 1];
      }
    }

    // We need __syncthreads here if the state is not aligned
    // because then we need to shift the vector [A, B, 0] -> [D, A, B]
    if (!aligned) { __syncthreads(); }

    // Store element in state
    states[idx] = prepared_element;
  }

#if CURVE_ID != 5
  template <typename S>
  __device__ __forceinline__ S sbox_alpha_five(S element)
  {
    S result = S::sqr(element);
    result = S::sqr(result);
    return result * element;
  }

  template <typename S, int T>
  __device__ S vecs_mul_matrix(S element, S* matrix, int element_number, int vec_number, S* shared_states)
  {
    __syncthreads();
    shared_states[threadIdx.x] = element;
    __syncthreads();

    typename S::Wide element_wide = S::mul_wide(shared_states[vec_number * T], matrix[element_number]);
#pragma unroll
    for (int i = 1; i < T; i++) {
      element_wide = element_wide + S::mul_wide(shared_states[vec_number * T + i], matrix[i * T + element_number]);
    }

    return S::reduce(element_wide);
  }

  template <typename S, int T>
  __device__ S full_round(
    S element,
    size_t rc_offset,
    int local_state_number,
    int element_number,
    bool multiply_by_mds,
    bool add_pre_round_constants,
    bool skip_rc,
    S* shared_states,
    const PoseidonConstants<S>& constants)
  {
    if (add_pre_round_constants) {
      element = element + constants.round_constants[rc_offset + element_number];
      rc_offset += T;
    }
    element = sbox_alpha_five(element);
    if (!skip_rc) { element = element + constants.round_constants[rc_offset + element_number]; }

    // Multiply all the states by mds matrix
    S* matrix = multiply_by_mds ? constants.mds_matrix : constants.non_sparse_matrix;
    return vecs_mul_matrix<S, T>(element, matrix, element_number, local_state_number, shared_states);
  }

  template <typename S, int T>
  __global__ void full_rounds(
    S* states, size_t number_of_states, size_t rc_offset, bool first_half, const PoseidonConstants<S> constants)
  {
    extern __shared__ S shared_states[];

    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    int state_number = idx / T;
    if (state_number >= number_of_states) { return; }
    int local_state_number = threadIdx.x / T;
    int element_number = idx % T;

    S new_el = states[idx];
    bool add_pre_round_constants = first_half;
    for (int i = 0; i < constants.full_rounds_half; i++) {
      new_el = full_round<S, T>(
        new_el, rc_offset, local_state_number, element_number, !first_half || (i < (constants.full_rounds_half - 1)),
        add_pre_round_constants, !first_half && (i == constants.full_rounds_half - 1), shared_states, constants);
      rc_offset += T;

      if (add_pre_round_constants) {
        rc_offset += T;
        add_pre_round_constants = false;
      }
    }
    states[idx] = new_el;
  }

  template <typename S, int T>
  __device__ S partial_round(S state[T], size_t rc_offset, int round_number, const PoseidonConstants<S>& constants)
  {
    S element = state[0];
    element = sbox_alpha_five(element);
    element = element + constants.round_constants[rc_offset];

    S* sparse_matrix = &constants.sparse_matrices[(T * 2 - 1) * round_number];

    typename S::Wide state_0_wide = S::mul_wide(element, sparse_matrix[0]);

#pragma unroll
    for (int i = 1; i < T; i++) {
      state_0_wide = state_0_wide + S::mul_wide(state[i], sparse_matrix[i]);
    }

    state[0] = S::reduce(state_0_wide);

#pragma unroll
    for (int i = 1; i < T; i++) {
      state[i] = state[i] + (element * sparse_matrix[T + i - 1]);
    }
  }

  template <typename S, int T>
  __global__ void
  partial_rounds(S* states, size_t number_of_states, size_t rc_offset, const PoseidonConstants<S> constants)
  {
    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (idx >= number_of_states) { return; }

    S state[T];
#pragma unroll
    for (int i = 0; i < T; i++) {
      state[i] = states[idx * T + i];
    }

    for (int i = 0; i < constants.partial_rounds; i++) {
      partial_round<S, T>(state, rc_offset, i, constants);
      rc_offset++;
    }

#pragma unroll
    for (int i = 0; i < T; i++) {
      states[idx * T + i] = state[i];
    }
  }

  // These function is just doing copy from the states to the output
  template <typename S, int T>
  __global__ void get_hash_results(S* states, size_t number_of_states, S* out)
  {
    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (idx >= number_of_states) { return; }

    out[idx] = states[idx * T + 1];
  }

  template <typename S, int T>
  __global__ void copy_recursive(S* state, size_t number_of_states, S* out)
  {
    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (idx >= number_of_states) { return; }

    state[(idx / (T - 1) * T) + (idx % (T - 1)) + 1] = out[idx];
  }
#else
  template <typename S>
  __device__ __forceinline__ S sbox_alpha_seven(S x)
  {
    S x2 = S::sqr(x);
    S x4 = S::sqr(x2);
    S x3 = x2 * x;
    return x3 * x4;
  }

  template <typename S, int T>
  __device__ S
  mds_row_shf(S element, S* mds_matrix_circ, S* mds_matrix_diag, int element_number, int vec_number, S* shared_states)
  {
    // S::print_debug("element_wide a = ", element, 2, "\n");
    // S::print_debug("element_wide b = ", mds_matrix_diag[element_number], 2, "\n");
    typename S::Wide element_wide{}; // = S::mul_wide(element, mds_matrix_diag[element_number]);

    __syncthreads();
    shared_states[threadIdx.x] = element;
    __syncthreads();

#pragma unroll
    for (int i = 0; i < T; i++) {
      // S::print_debug("element_wide a = ", shared_states[vec_number * T + (i + element_number) % T], 2, "\n");
      // S::print_debug("element_wide b = ", mds_matrix_circ[i], 2, "\n");
      element_wide =
        element_wide + S::mul_wide(shared_states[vec_number * T + (i + element_number) % T], mds_matrix_circ[i]);
      // S::print_debug("element_wide = ", element_wide, 4, "\n");
    }
    element_wide = element_wide + S::mul_wide(element, mds_matrix_diag[element_number]);
    // S::print_debug("element_wide before reduce = ", element_wide, 4, "\n");
    S res = S::reduce96(element_wide);
    // S::print_debug("element_wide after reduce = ", res, 2, "\n");
    return res;
  }

  template <typename S, int T>
  __device__ void print_states_device(S* states)
  {
    printf("states: ");
    for (int i = 0; i < T; i++) {
      S::print_debug("", states[i], 2, ", ");
    }
    printf("\n");
  }

  template <typename S, int T>
  __global__ void print_states(S* states, size_t number_of_states)
  {
    for (int idx = 0; idx < number_of_states; idx++) {
      printf("states %d : ", idx);
      for (int i = 0; i < T; i++) {
        S::print_debug("", states[idx * T + i], 2, ", ");
      }
      printf("\n");
    }
  }

  template <typename S, int T>
  __device__ S full_round(
    S element,
    size_t rc_offset,
    int local_state_number,
    int element_number,
    S* shared_states,
    const PoseidonConstants<S>& constants)
  {
    // S::print_debug("element = ", element, 2, "\n");
    element = element + constants.all_round_constants[rc_offset * T + element_number];
    element = sbox_alpha_seven(element);
    // S::print_debug("element = ", element, 2, "\n");

    // Multiply all the states by mds matrix
    S* mds_matrix_circ = constants.mds_matrix_circ;
    S* mds_matrix_diag = constants.mds_matrix_diag;
    S result =
      mds_row_shf<S, T>(element, mds_matrix_circ, mds_matrix_diag, element_number, local_state_number, shared_states);
    // S::print_debug("result = ", result, 2, "\n");
    return result;
  }

  template <typename S, int T>
  __global__ void
  full_rounds(S* states, size_t number_of_states, size_t rc_offset, const PoseidonConstants<S> constants)
  {
    extern __shared__ S shared_states[];

    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    int state_number = idx / T;
    if (state_number >= number_of_states) { return; }
    // printf("state_number %d\n", state_number);
    // printf("T %d\n", T);
    // printf("number_of_states %d\n", number_of_states);
    // printf("full_rounds\n");
    int local_state_number = threadIdx.x / T;
    int element_number = idx % T;

    // printf("full_rounds 0\n");
    S new_el = states[idx];
    // printf("full_rounds 1\n");
    for (int i = 0; i < constants.full_rounds_half; i++) {
      new_el = full_round<S, T>(new_el, rc_offset, local_state_number, element_number, shared_states, constants);
      // printf("full_rounds 2\n");
      rc_offset += 1;
    }
    // printf("full_rounds 3\n");
    states[idx] = new_el;
    // printf("full_rounds done\n");
  }

  template <typename S, int T>
  __device__ S partial_round(S state[T], int round_number, const PoseidonConstants<S>& constants)
  {
    S element = state[0];
    element = sbox_alpha_seven(element);
    element = element + constants.fast_partial_round_constants[round_number];

    typename S::U160 d_sum{};

#pragma unroll
    for (int i = 1; i < T; i++) {
      S t = constants.fast_partial_round_w_hats[round_number * constants.fast_partial_round_w_hats_len_y + i - 1];
      S::add_limbs_device(d_sum.limbs_storage, S::mul_wide(t, state[i]).limbs_storage, d_sum.limbs_storage);
    }
    S mds0to0 = S::from(constants.mds0to0);
    S::add_limbs_device(d_sum.limbs_storage, S::mul_wide(mds0to0, element).limbs_storage, d_sum.limbs_storage);

#pragma unroll
    for (int i = 1; i < T; i++) {
      state[i] =
        state[i] +
        (element * constants.fast_partial_round_vs[(round_number * constants.fast_partial_round_vs_len_y) + (i - 1)]);
    }
    // S::print_debug("d_sum = ", d_sum, 5, "\n");
    state[0] = S::reduce160(d_sum);
  }

  template <typename S, int T>
  __global__ void partial_rounds(S* states, size_t number_of_states, const PoseidonConstants<S> constants)
  {
    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (idx >= number_of_states) { return; }

    S state[T];
#pragma unroll
    for (int i = 0; i < T; i++) {
      state[i] = states[idx * T + i] + constants.fast_partial_first_round_constant[i];
    }

    S state_mds[T] = {};
    state_mds[0] = state[0];

#pragma unroll
    for (int r = 1; r < T; r++) {
      for (int c = 1; c < T; c++) {
        S t = constants
                .fast_partial_round_initial_matrix[(r - 1) * constants.fast_partial_round_initial_matrix_len_y + c - 1];
        state_mds[c] = state_mds[c] + state[r] * t;
      }
    }

    for (int i = 0; i < constants.partial_rounds; i++) {
      partial_round<S, T>(state_mds, i, constants);
    }

#pragma unroll
    for (int i = 0; i < T; i++) {
      states[idx * T + i] = state_mds[i];
    }
  }

  // These function is just doing copy from the states to the output
  template <typename S, int T>
  __global__ void get_hash_results(S* states, size_t number_of_states, S* out)
  {
    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (idx >= number_of_states) { return; }

    for (int i = 0; i < 4; i++) {
      out[idx * 4 + i] = states[idx * T + i];
    }
  }

  template <typename S, int T>
  __global__ void copy_recursive(S* state, size_t number_of_states, S* out)
  {
    int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
    if (idx >= number_of_states) { return; }

    for (int j = 0; j < 12; j++) {
      if (j < 8)
        state[idx * 12 + j] = out[idx * 8 + j];
      else
        state[idx * 12 + j] = S::zero();
    }
  }
#endif
} // namespace poseidon