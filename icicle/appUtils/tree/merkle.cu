#include "merkle.cuh"
#define CURVE_ID 5
namespace merkle {
#if CURVE_ID != 5
  /// Flattens the tree digests and sum them up to get
  /// the memory needed to contain all the digests
  template <typename S>
  size_t get_digests_len(uint32_t height, uint32_t arity)
  {
    size_t digests_len = 0;
    size_t row_length = 1;
    for (int i = 1; i < height; i++) {
      digests_len += row_length;
      row_length *= arity;
    }

    return digests_len;
  }

  /// Constructs merkle subtree without parallelization
  /// The digests are aligned sequentially per row
  /// Example:
  ///
  /// Big tree:
  ///
  ///        1
  ///       / \
  ///      2   3
  ///     / \ / \
  ///    4  5 6  7
  ///
  /// Subtree 1    Subtree 2
  ///    2            3
  ///   / \          / \
  ///  4   5        6   7
  ///
  /// Digests array for subtree 1:
  /// [4 5 . . 2 . .]
  /// |   |    |
  /// -----    V
  ///   |    Segment (offset = 4, subtree_idx = 0)
  ///   v
  /// Segment (offset = 0, subtree_idx = 0)
  ///
  /// Digests array for subtree 2:
  /// [. . 6 7 . 3 .]
  ///     |   |
  ///     -----
  ///       |
  ///       v
  ///    Segment (offset = 0, subtree_idx = 1)
  ///
  /// Total digests array:
  /// [4 5 6 7 2 3 .]
  template <typename S, int T>
  cudaError_t build_merkle_subtree(
    S* state,
    S* digests,
    size_t subtree_idx,
    size_t subtree_height,
    S* big_tree_digests,
    size_t start_segment_size,
    size_t start_segment_offset,
    int keep_rows,
    const PoseidonConstants<S>& poseidon,
    cudaStream_t& stream)
  {
    int arity = T - 1;

    PoseidonConfig config = default_poseidon_config<S>(T);
    config.are_inputs_on_device = true;
    config.are_outputs_on_device = true;
    config.input_is_a_state = true;
    config.loop_state = true;
    config.ctx.stream = stream;

    size_t leaves_size = pow(arity, subtree_height - 1);
    uint32_t number_of_blocks = leaves_size / arity;
    size_t segment_size = start_segment_size;
    size_t segment_offset = start_segment_offset;

    while (number_of_blocks > 0) {
      cudaError_t poseidon_res = poseidon_hash<S, T>(state, digests, number_of_blocks, poseidon, config);
      CHK_IF_RETURN(poseidon_res);

      if (!keep_rows || subtree_height <= keep_rows + 1) {
        S* digests_with_offset = big_tree_digests + segment_offset + subtree_idx * number_of_blocks;
        CHK_IF_RETURN(
          cudaMemcpyAsync(digests_with_offset, digests, number_of_blocks * sizeof(S), cudaMemcpyDeviceToHost, stream));
        segment_offset += segment_size;
      }

      segment_size /= arity;
      subtree_height--;
      number_of_blocks /= arity;
      config.aligned = true;
    }

    return CHK_LAST();
  }

  template <typename S, int T>
  cudaError_t build_merkle_tree(
    const S* leaves,
    S* digests,
    uint32_t height,
    const poseidon::PoseidonConstants<S>& poseidon,
    const TreeBuilderConfig& config)
  {
    CHK_INIT_IF_RETURN();
    cudaStream_t& stream = config.ctx.stream;

    int arity = T - 1;
    uint32_t number_of_leaves = pow(arity, (height - 1));

    // This will determine how much splitting do we need to do
    // `number_of_streams` subtrees should fit in the device
    // This means each subtree should fit in `STREAM_CHUNK_SIZE` memory
    uint32_t number_of_subtrees = 1;
    uint32_t subtree_height = height;
    uint32_t subtree_leaves_size = pow(arity, height - 1);
    uint32_t subtree_state_size = subtree_leaves_size / arity * T;
    uint32_t subtree_digests_size = get_digests_len<S>(subtree_height, arity);
    size_t subtree_memory_required = sizeof(S) * (subtree_state_size + subtree_digests_size);
    while (subtree_memory_required > STREAM_CHUNK_SIZE) {
      number_of_subtrees *= arity;
      subtree_height--;
      subtree_leaves_size /= arity;
      subtree_state_size = subtree_leaves_size / arity * T;
      subtree_digests_size = subtree_state_size / arity;
      subtree_memory_required = sizeof(S) * (subtree_state_size + subtree_digests_size);
    }
    int cap_height = height - subtree_height + 1;
    size_t caps_len = pow(arity, cap_height - 1);

    size_t available_memory, _total_memory;
    CHK_IF_RETURN(cudaMemGetInfo(&available_memory, &_total_memory));
    available_memory -= GIGA / 8; // Leave 128 MB

    // We can effectively parallelize memory copy with streams
    // as long as they don't operate on more than `STREAM_CHUNK_SIZE` bytes
    const size_t number_of_streams = std::min((uint32_t)(available_memory / STREAM_CHUNK_SIZE), number_of_subtrees);
    cudaStream_t* streams = static_cast<cudaStream_t*>(malloc(sizeof(cudaStream_t) * number_of_streams));
    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamCreate(&streams[i]));
    }

#if !defined(__CUDA_ARCH__) && defined(MERKLE_DEBUG)
    std::cout << "Available memory = " << available_memory / 1024 / 1024 << " MB" << std::endl;
    std::cout << "Number of streams = " << number_of_streams << std::endl;
    std::cout << "Number of subtrees = " << number_of_subtrees << std::endl;
    std::cout << "Height of a subtree = " << subtree_height << std::endl;
    std::cout << "Cutoff height = " << height - subtree_height + 1 << std::endl;
    std::cout << "Number of leaves in a subtree = " << subtree_leaves_size << std::endl;
    std::cout << "State of a subtree = " << subtree_state_size << std::endl;
    std::cout << "Digest elements for a subtree = " << get_digests_len<S>(subtree_height, arity) << std::endl;
    std::cout << "Size of 1 subtree states = " << subtree_state_size * sizeof(S) / 1024 / 1024 << " MB" << std::endl;
    std::cout << "Size of 1 subtree digests = " << subtree_digests_size * sizeof(S) / 1024 / 1024 << " MB" << std::endl;
#endif

    // Allocate memory for the leaves and digests
    // These are shared by streams in a pool
    S *states_ptr, *digests_ptr;
    CHK_IF_RETURN(cudaMallocAsync(&states_ptr, subtree_state_size * number_of_streams * sizeof(S), stream))
    CHK_IF_RETURN(cudaMallocAsync(&digests_ptr, subtree_digests_size * number_of_streams * sizeof(S), stream))
    // Wait for these allocations to finish
    CHK_IF_RETURN(cudaStreamSynchronize(stream));

    bool caps_mode = config.keep_rows && config.keep_rows < cap_height;
    S* caps;
    if (caps_mode) { caps = static_cast<S*>(malloc(caps_len * sizeof(S))); }

    for (size_t subtree_idx = 0; subtree_idx < number_of_subtrees; subtree_idx++) {
      size_t stream_idx = subtree_idx % number_of_streams;
      cudaStream_t subtree_stream = streams[stream_idx];

      const S* subtree_leaves = leaves + subtree_idx * subtree_leaves_size;
      S* subtree_state = states_ptr + stream_idx * subtree_state_size;
      S* subtree_digests = digests_ptr + stream_idx * subtree_digests_size;

      // Copy the first level from RAM / device to device
      // The pitch property of cudaMemcpy2D resolves shape differences
      CHK_IF_RETURN(cudaMemcpy2DAsync(
        subtree_state, T * sizeof(S),      // Device pointer and device pitch
        subtree_leaves, arity * sizeof(S), // Host pointer and pitch
        arity * sizeof(S),                 // Size of the source matrix (Arity)
        subtree_leaves_size / arity,       // Size of the source matrix (Number of blocks)
        config.are_inputs_on_device ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice, subtree_stream));

      int subtree_keep_rows = 0;
      if (config.keep_rows) {
        int diff = config.keep_rows - cap_height + 1;
        subtree_keep_rows = diff <= 0 ? 1 : diff;
      }
      size_t start_segment_size = number_of_leaves / arity;
      cudaError_t subtree_result = build_merkle_subtree<S, T>(
        subtree_state,              // state
        subtree_digests,            // digests
        subtree_idx,                // subtree_idx
        subtree_height,             // subtree_height
        caps_mode ? caps : digests, // big_tree_digests
        start_segment_size,         // start_segment_size
        0,                          // start_segment_offset
        subtree_keep_rows,          // keep_rows
        poseidon,                   // hash
        subtree_stream              // stream
      );
      CHK_IF_RETURN(subtree_result);
    }

    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamSynchronize(streams[i]));
    }

    // Finish the top-level tree if any
    if (cap_height > 1) {
      size_t start_segment_size = caps_len / arity;
      size_t start_segment_offset = 0;
      if (!caps_mode) {
        size_t layer_size = pow(arity, config.keep_rows - 1);
        for (int i = 0; i < config.keep_rows - cap_height + 1; i++) {
          start_segment_offset += layer_size;
          layer_size /= arity;
        }
      }
      CHK_IF_RETURN(cudaMemcpy2DAsync(
        states_ptr, T * sizeof(S), caps_mode ? caps : (digests + start_segment_offset - caps_len), arity * sizeof(S),
        arity * sizeof(S),
        caps_len / arity,                 // Size of the source
        cudaMemcpyHostToDevice, stream)); // Direction and stream

      cudaError_t top_tree_result = build_merkle_subtree<S, T>(
        states_ptr,           // state
        digests_ptr,          // digests
        0,                    // subtree_idx
        cap_height,           // subtree_height
        digests,              // big_tree_digests
        start_segment_size,   // start_segment_size
        start_segment_offset, // start_segment_offset
        config.keep_rows,     // keep_rows
        poseidon,             // hash
        stream                // stream
      );
      CHK_IF_RETURN(top_tree_result);
      if (caps_mode) { free(caps); }
    }

    CHK_IF_RETURN(cudaFreeAsync(states_ptr, stream));
    CHK_IF_RETURN(cudaFreeAsync(digests_ptr, stream));
    if (!config.is_async) return CHK_STICKY(cudaStreamSynchronize(stream));
    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamSynchronize(streams[i]));
      CHK_IF_RETURN(cudaStreamDestroy(streams[i]));
    }
    free(streams);
    return CHK_LAST();
  }

  extern "C" cudaError_t CONCAT_EXPAND(CURVE, BuildPoseidonMerkleTree)(
    const curve_config::scalar_t* leaves,
    curve_config::scalar_t* digests,
    uint32_t height,
    int arity,
    PoseidonConstants<curve_config::scalar_t>& constants,
    TreeBuilderConfig& config)
  {
    switch (arity) {
    case 2:
      return build_merkle_tree<curve_config::scalar_t, 3>(leaves, digests, height, constants, config);
    case 4:
      return build_merkle_tree<curve_config::scalar_t, 5>(leaves, digests, height, constants, config);
    case 8:
      return build_merkle_tree<curve_config::scalar_t, 9>(leaves, digests, height, constants, config);
    case 11:
      return build_merkle_tree<curve_config::scalar_t, 12>(leaves, digests, height, constants, config);
    default:
      THROW_ICICLE_ERR(IcicleError_t::InvalidArgument, "BuildPoseidonMerkleTree: #arity must be one of [2, 4, 8, 11]");
    }
    return CHK_LAST();
  }

#else
#define RATE 8
  // transform the addr format to the plonky2 tree format

  size_t addr_to_index_fast(uint32_t level, uint32_t idx, uint32_t height, uint32_t cap_height)
  {
    uint32_t subtree_height = height - cap_height;
    if (subtree_height == 1) return idx;

    uint32_t number_of_leaves = 1 << (height - 1);
    uint32_t number_of_nodes = 1 << (height - level - 1);
    uint32_t number_of_nodes_in_subtree = number_of_nodes >> cap_height;
    uint32_t subtree_idx = idx / number_of_nodes_in_subtree;

    uint32_t num_digests = 2 * (number_of_leaves - (1 << cap_height));
    uint32_t digest_buf_offset = (num_digests >> cap_height) * subtree_idx;

    uint32_t number_of_left_nodes = (idx & 1) ? (1 << (level + 1)) - 1 : (1 << (level + 1)) - 2;
    uint32_t mask = (1 << (subtree_height - 1)) - 1;
    number_of_left_nodes += (((idx >> 1) << (level + 1)) & mask) << 1;

    return digest_buf_offset + number_of_left_nodes;
  }

  template <typename S, int T>
  cudaError_t build_merkle_subtree(
    S* state,
    S* digests,
    size_t subtree_idx,
    size_t height,
    size_t cap_height,
    S* big_tree_digests,
    size_t start_segment_size,
    size_t start_segment_offset,
    size_t digest_length,
    const PoseidonConstants<S>& poseidon,
    cudaStream_t& stream)
  {
    PoseidonConfig config = default_poseidon_config<S>(T);
    config.are_inputs_on_device = true;
    config.are_outputs_on_device = true;
    config.input_is_a_state = true;
    config.loop_state = true;
    config.ctx.stream = stream;

    size_t subtree_height = height - cap_height;
    size_t leaves_size = pow(2, subtree_height - 1);
    size_t number_of_blocks = leaves_size / 2;
    size_t segment_size = start_segment_size;
    size_t segment_offset = start_segment_offset;
    size_t number_of_leaves = 1 << (height - 1);
    size_t num_digests = 2 * (number_of_leaves - (1 << cap_height));

    size_t level = 1;

    while (number_of_blocks > 0) {
      if (number_of_blocks == 1) config.loop_state = false;
      cudaError_t poseidon_res = poseidon_hash<S, T>(state, digests, number_of_blocks, poseidon, config);
      CHK_IF_RETURN(poseidon_res);

      for (size_t i = 0; i < number_of_blocks; i++) {
        size_t addr;
        if (number_of_blocks > 1)
          addr = addr_to_index_fast(level, subtree_idx * number_of_blocks + i, height, cap_height);
        else
          addr = num_digests + subtree_idx;

        S* digests_with_offset = big_tree_digests + addr * digest_length;
        CHK_IF_RETURN(cudaMemcpyAsync(
          digests_with_offset, digests + i * digest_length, digest_length * sizeof(S), cudaMemcpyDeviceToHost, stream))
      }

      segment_offset += segment_size;
      segment_size /= 2;
      subtree_height--;
      number_of_blocks /= 2;
      level++;
    }

    return CHK_LAST();
  }

  template <typename S, int T>
  cudaError_t absorb_leaves(
    const S* leaves,
    S* digests,
    uint32_t height,
    uint32_t cap_height,
    uint32_t leaf_length,
    uint32_t digest_length,
    const poseidon::PoseidonConstants<S>& poseidon,
    const TreeBuilderConfig& config)
  {
    cudaStream_t& stream = config.ctx.stream;

    uint32_t number_of_leaves = pow(2, (height - 1));

    uint32_t chunk_size =
      std::min(uint32_t((STREAM_CHUNK_SIZE + T * sizeof(S) - 1) / (T * sizeof(S))), number_of_leaves);
    uint32_t num_chunks = (number_of_leaves + chunk_size - 1) / chunk_size;

    size_t available_memory, _total_memory;
    CHK_IF_RETURN(cudaMemGetInfo(&available_memory, &_total_memory));
    available_memory -= GIGA / 8; // Leave 128 MB

    const size_t number_of_streams = std::min((uint32_t)(available_memory / STREAM_CHUNK_SIZE), num_chunks);
    cudaStream_t* streams = static_cast<cudaStream_t*>(malloc(sizeof(cudaStream_t) * number_of_streams));
    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamCreate(&streams[i]));
    }

    S *states_ptr, *digests_ptr;
    CHK_IF_RETURN(cudaMallocAsync(&states_ptr, chunk_size * T * number_of_streams * sizeof(S), stream))
    CHK_IF_RETURN(cudaMemsetAsync(states_ptr, 0, chunk_size * T * number_of_streams * sizeof(S), stream))
    CHK_IF_RETURN(cudaMallocAsync(&digests_ptr, chunk_size * digest_length * number_of_streams * sizeof(S), stream))
    // Wait for these allocations to finish
    CHK_IF_RETURN(cudaStreamSynchronize(stream));

#if !defined(__CUDA_ARCH__) && defined(MERKLE_DEBUG)
    std::cout << "Available memory = " << available_memory / 1024 / 1024 << " MB" << std::endl;
    std::cout << "Number of streams = " << number_of_streams << std::endl;
    std::cout << "Chunk size = " << chunk_size << std::endl;
    std::cout << "Number of chunks = " << num_chunks << std::endl;
#endif

    for (size_t leaf_idx = 0; leaf_idx < leaf_length; leaf_idx += RATE) {
      for (size_t chunk_idx = 0; chunk_idx < num_chunks; chunk_idx++) {
        size_t stream_idx = chunk_idx % number_of_streams;
        cudaStream_t sub_stream = streams[stream_idx];
        size_t copy_length = std::min((uint32_t)RATE, uint32_t(leaf_length - leaf_idx));
        size_t chunk_size_actual = std::min((uint32_t)chunk_size, uint32_t(number_of_leaves - chunk_idx * chunk_size));

        const S* leaves_chunk_ptr = leaves + chunk_idx * chunk_size * leaf_length;
        S* states_chunk_ptr = states_ptr + stream_idx * chunk_size * T;
        S* digests_chunk_ptr = digests_ptr + stream_idx * chunk_size * digest_length;

        // CHK_IF_RETURN(cudaMemcpy2DAsync(
        //   states_chunk_ptr, T * sizeof(S),      // Device pointer and device pitch
        //   leaves_chunk_ptr + leaf_idx, leaf_length * sizeof(S), // Host pointer and pitch
        //   copy_length * sizeof(S),                 // Size of the source matrix
        //   chunk_size_actual,       // Size of the source matrix
        //   config.are_inputs_on_device ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice, subtree_stream));
        for (size_t i = 0; i < chunk_size_actual; i++) {
          CHK_IF_RETURN(cudaMemcpyAsync(
            states_chunk_ptr + i * T, leaves_chunk_ptr + i * leaf_length + leaf_idx, copy_length * sizeof(S),
            config.are_inputs_on_device ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice, sub_stream));
        }

        PoseidonConfig pconfig = default_poseidon_config<S>(T);
        pconfig.are_inputs_on_device = true;
        pconfig.are_outputs_on_device = true;
        pconfig.input_is_a_state = true;
        pconfig.loop_state = false;
        pconfig.ctx.stream = sub_stream;

        cudaError_t poseidon_res =
          poseidon_hash<S, T>(states_chunk_ptr, digests_chunk_ptr, chunk_size_actual, poseidon, pconfig);
        CHK_IF_RETURN(poseidon_res);

        if (leaf_idx + RATE >= leaf_length) {
          for (size_t i = 0; i < chunk_size_actual; i++) {
            size_t addr = addr_to_index_fast(0, chunk_idx * chunk_size + i, height, cap_height);
            S* digests_with_offset = digests + addr * digest_length;
            CHK_IF_RETURN(cudaMemcpyAsync(
              digests_with_offset, digests_chunk_ptr + i * digest_length, digest_length * sizeof(S),
              cudaMemcpyDeviceToHost, sub_stream))
          }
        }
      }
    }

    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamSynchronize(streams[i]));
    }

    CHK_IF_RETURN(cudaFreeAsync(states_ptr, stream));
    CHK_IF_RETURN(cudaFreeAsync(digests_ptr, stream));
    if (!config.is_async) return CHK_STICKY(cudaStreamSynchronize(stream));
    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamSynchronize(streams[i]));
      CHK_IF_RETURN(cudaStreamDestroy(streams[i]));
    }
    free(streams);
    return CHK_LAST();
  }

  template <typename S, int T>
  cudaError_t split_merkle_tree(
    S* digests,
    uint32_t height,
    uint32_t cap_height,
    uint32_t digest_length,
    const poseidon::PoseidonConstants<S>& poseidon,
    const TreeBuilderConfig& config)
  {
    cudaStream_t& stream = config.ctx.stream;

    uint32_t number_of_leaves = pow(2, (height - 1));
    uint32_t number_of_subtrees = 1 << cap_height;
    uint32_t subtree_height = height - cap_height;
    uint32_t subtree_leaves_size = pow(2, subtree_height - 1);
    uint32_t subtree_state_size = subtree_leaves_size / 2 * T;
    uint32_t subtree_digests_size = ((2 * number_of_leaves - (1 << cap_height)) >> cap_height) * digest_length;
    size_t subtree_memory_required = sizeof(S) * (subtree_state_size + subtree_digests_size);
    assert(subtree_memory_required <= STREAM_CHUNK_SIZE);

    size_t available_memory, _total_memory;
    CHK_IF_RETURN(cudaMemGetInfo(&available_memory, &_total_memory));
    available_memory -= GIGA / 8; // Leave 128 MB

    const size_t number_of_streams = std::min((uint32_t)(available_memory / STREAM_CHUNK_SIZE), number_of_subtrees);
    cudaStream_t* streams = static_cast<cudaStream_t*>(malloc(sizeof(cudaStream_t) * number_of_streams));
    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamCreate(&streams[i]));
    }

#if !defined(__CUDA_ARCH__) && defined(MERKLE_DEBUG)
    std::cout << "Available memory = " << available_memory / 1024 / 1024 << " MB" << std::endl;
    std::cout << "Number of streams = " << number_of_streams << std::endl;
    std::cout << "Number of subtrees = " << number_of_subtrees << std::endl;
    std::cout << "Height of a subtree = " << subtree_height << std::endl;
    std::cout << "Number of leaves in a subtree = " << subtree_leaves_size << std::endl;
    std::cout << "State of a subtree = " << subtree_state_size << std::endl;
    std::cout << "Digest elements for a subtree = " << ((2 * (number_of_leaves - (1 << cap_height)) >> cap_height) + 1)
              << std::endl;
    std::cout << "Size of 1 subtree states = " << subtree_state_size * sizeof(S) / 1024 / 1024 << " MB" << std::endl;
    std::cout << "Size of 1 subtree digests = " << subtree_digests_size * sizeof(S) / 1024 / 1024 << " MB" << std::endl;
#endif

    // Allocate memory for the leaves and digests
    // These are shared by streams in a pool
    S *states_ptr, *digests_ptr;
    CHK_IF_RETURN(cudaMallocAsync(&states_ptr, subtree_state_size * number_of_streams * sizeof(S), stream))
    CHK_IF_RETURN(cudaMemsetAsync(states_ptr, 0, subtree_state_size * number_of_streams * sizeof(S), stream))
    CHK_IF_RETURN(cudaMallocAsync(&digests_ptr, subtree_digests_size * number_of_streams * sizeof(S), stream))
    // Wait for these allocations to finish
    CHK_IF_RETURN(cudaStreamSynchronize(stream));

    for (size_t subtree_idx = 0; subtree_idx < number_of_subtrees; subtree_idx++) {
      size_t stream_idx = subtree_idx % number_of_streams;
      cudaStream_t subtree_stream = streams[stream_idx];

      S* subtree_state = states_ptr + stream_idx * subtree_state_size;
      S* subtree_digests = digests_ptr + stream_idx * subtree_digests_size;

      // Copy the first level from RAM / device to device
      // The pitch property of cudaMemcpy2D resolves shape differences
      for (size_t i = 0; i < pow(2, subtree_height - 1); i++) {
        uint32_t idx = i + (subtree_idx << (subtree_height - 1));
        size_t addr = addr_to_index_fast(0, idx, height, cap_height);
        CHK_IF_RETURN(cudaMemcpyAsync(
          subtree_state + i / 2 * T + i % 2 * digest_length, digests + addr * digest_length, digest_length * sizeof(S),
          config.are_inputs_on_device ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice, subtree_stream));
      }
      // CHK_IF_RETURN(cudaMemcpy2DAsync(
      //   subtree_state, T * sizeof(S),                  // Device pointer and device pitch
      //   subtree_leaves, 2 * digest_length * sizeof(S), // Host pointer and pitch
      //   2 * digest_length * sizeof(S),                 // Size of the source matrix (Arity)
      //   subtree_leaves_size,                           // Size of the source matrix (Number of blocks)
      //   config.are_inputs_on_device ? cudaMemcpyDeviceToDevice : cudaMemcpyHostToDevice, subtree_stream));

      size_t start_segment_size = number_of_leaves;
      cudaError_t subtree_result = build_merkle_subtree<S, T>(
        subtree_state,      // state
        subtree_digests,    // digests
        subtree_idx,        // subtree_idx
        height,             // height
        cap_height,         // cap_height
        digests,            // big_tree_digests
        start_segment_size, // start_segment_size
        0,                  // start_segment_offset
        digest_length,      // digest_length
        poseidon,           // hash
        subtree_stream      // stream
      );
      CHK_IF_RETURN(subtree_result);
    }

    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamSynchronize(streams[i]));
    }

    CHK_IF_RETURN(cudaFreeAsync(states_ptr, stream));
    CHK_IF_RETURN(cudaFreeAsync(digests_ptr, stream));
    if (!config.is_async) return CHK_STICKY(cudaStreamSynchronize(stream));
    for (size_t i = 0; i < number_of_streams; i++) {
      CHK_IF_RETURN(cudaStreamSynchronize(streams[i]));
      CHK_IF_RETURN(cudaStreamDestroy(streams[i]));
    }
    free(streams);
    return CHK_LAST();
  }

  template <typename S, int T>
  cudaError_t build_merkle_tree(
    const S* leaves,
    S* digests,
    uint32_t height,
    uint32_t cap_height,
    uint32_t leaf_length,
    uint32_t digest_length,
    const poseidon::PoseidonConstants<S>& poseidon,
    const TreeBuilderConfig& config)
  {
    CHK_INIT_IF_RETURN();
    cudaError_t absorb_res =
      absorb_leaves<S, T>(leaves, digests, height, cap_height, leaf_length, digest_length, poseidon, config);
    CHK_IF_RETURN(absorb_res);
    uint32_t number_of_leaves = pow(2, (height - 1));
    if (cap_height < height - 1) {
      // The leaves are actually the digests in the first level.
      cudaError_t subtree_res = split_merkle_tree<S, T>(digests, height, cap_height, digest_length, poseidon, config);
      CHK_IF_RETURN(subtree_res);
    }
    return CHK_LAST();
  }
#endif
} // namespace merkle