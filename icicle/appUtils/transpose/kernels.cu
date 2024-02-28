#include "transpose.cuh"

namespace transpose {
#define INDEX(ROW, COL, INNER) ((ROW) * (INNER) + (COL))

  // template <typename S, int T>
  // __global__ void
  // full_rounds(S* states, size_t number_of_states, size_t rc_offset, const PoseidonConstants<S> constants)
  // {
  //   extern __shared__ S shared_states[];

  //   int idx = (blockIdx.x * blockDim.x) + threadIdx.x;
  //   int state_number = idx / T;
  //   if (state_number >= number_of_states) { return; }
  //   // printf("state_number %d\n", state_number);
  //   // printf("T %d\n", T);
  //   // printf("number_of_states %d\n", number_of_states);
  //   // printf("full_rounds\n");
  //   int local_state_number = threadIdx.x / T;
  //   int element_number = idx % T;

  //   // printf("full_rounds 0\n");
  //   S new_el = states[idx];
  //   // printf("full_rounds 1\n");
  //   for (int i = 0; i < constants.full_rounds_half; i++) {
  //     new_el = full_round<S, T>(new_el, rc_offset, local_state_number, element_number, shared_states, constants);
  //     // printf("full_rounds 2\n");
  //     rc_offset += 1;
  //   }
  //   // printf("full_rounds 3\n");
  //   states[idx] = new_el;
  //   // printf("full_rounds done\n");
  // }
  template <typename S>
  __global__ void transposeSmemUnrollPadDyn(const S* src, S* dst, const int nrows, const int ncols)
  {
    // dynamic shared memory
    extern __shared__ S tile[];

    uint32_t row = blockIdx.y * blockDim.y + threadIdx.y;
    uint32_t col = (2 * blockIdx.x * blockDim.x) + threadIdx.x;

    uint32_t row2 = row;
    uint32_t col2 = col + blockDim.x;

    // linear global memory index for original matrix
    uint32_t offset = INDEX(row, col, ncols);
    uint32_t offset2 = INDEX(row2, col2, ncols);

    // thread index in transposed block
    uint32_t bidx = threadIdx.y * blockDim.x + threadIdx.x;
    uint32_t irow = bidx / blockDim.y;
    uint32_t icol = bidx % blockDim.y;

    // coordinate in transposed matrix
    uint32_t transposed_offset = INDEX(col, row, nrows);
    uint32_t transposed_offset2 = INDEX(col2, row2, nrows);
    // printf(
    //   "\
    //   nrows: %d, ncols: %d\n \
    //   threadIdx:(%d, %d, %d)\n \
    //   blockIdx:(%d, %d, %d)\n \
    //   row: %d, col: %d, offset: %d\n",
    //   nrows, ncols, threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y, blockIdx.z, row, col, offset);

    if (row < nrows && col < ncols) {
      // printf("row: %d, col: %d, offset: %d\n", row, col, offset);
      tile[INDEX(threadIdx.y, threadIdx.x, BDIMX * 2 + IPAD)] = src[offset];
    }
    if (row2 < nrows && col2 < ncols) {
      // printf("row2: %d, col2: %d, offset2: %d\n", row2, col2, offset2);
      tile[INDEX(threadIdx.y, blockDim.x + threadIdx.x, BDIMX * 2 + IPAD)] = src[offset2];
    }

    __syncthreads();

    if (row < nrows && col < ncols) {
      // printf("irow: %d, icol: %d, sm offset: %d\n", irow, icol, INDEX(irow, icol, BDIMX * 2 + IPAD));
      dst[transposed_offset] = tile[INDEX(irow, icol, BDIMX * 2 + IPAD)];
    }
    if (row2 < nrows && col2 < ncols) {
      // printf(
      // "irow: %d, icol: %d, sm offset: %d\n", irow, blockDim.x + icol,
      // INDEX(irow, blockDim.x + icol, BDIMX * 2 + IPAD));
      dst[transposed_offset2] = tile[INDEX(irow, blockDim.x + icol, BDIMX * 2 + IPAD)];
    }
  }

} // namespace transpose