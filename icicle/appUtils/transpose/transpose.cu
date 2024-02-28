#include "transpose.cuh"
#include "kernels.cu"

namespace transpose {

  template <typename S>
  cudaError_t transpose_mem(S* input, S* output, const TransposeConfig& config)
  {
    CHK_INIT_IF_RETURN();
    cudaStream_t& stream = config.ctx.stream;
    S* src;
    if (config.are_inputs_on_device) {
      src = input;
    } else {
      printf("matrix size in bytes: %lu\n", TKC<S>::size_bytes(config));
      CHK_IF_RETURN(cudaMallocAsync(&src, TKC<S>::size_bytes(config), stream))
      CHK_IF_RETURN(cudaMemcpyAsync(src, input, TKC<S>::size_bytes(config), cudaMemcpyHostToDevice, stream));
    }

    S* dst;
    if (config.are_outputs_on_device) {
      dst = output;
    } else {
      CHK_IF_RETURN(cudaMallocAsync(&dst, TKC<S>::size_bytes(config), stream))
    }

    dim3 grid = TKC<S>::get_grid(config);
    dim3 block = TKC<S>::get_block();

    // transposeSmemUnrollPadDynP<<<grid, block, 0, stream>>>();
    transposeSmemUnrollPadDyn<S><<<grid, block, TKC<S>::sm_size_bytes(), stream>>>(src, dst, config.nrows, config.ncols);
    // checkaa<8><<<grid, block>>>();

    if (!config.are_inputs_on_device) CHK_IF_RETURN(cudaFreeAsync(src, stream));

    if (!config.are_outputs_on_device) {
      CHK_IF_RETURN(cudaMemcpyAsync(output, dst, TKC<S>::size_bytes(config), cudaMemcpyDeviceToHost, stream));
      CHK_IF_RETURN(cudaFreeAsync(dst, stream));
    }

    if (!config.is_async) return CHK_STICKY(cudaStreamSynchronize(stream));
    return CHK_LAST();
  }

  extern "C" cudaError_t CONCAT_EXPAND(CURVE, TransposeMem)(
    curve_config::scalar_t* input, curve_config::scalar_t* output, TransposeConfig& config)
  {
    transpose_mem<curve_config::scalar_t>(input, output, config);
    return CHK_LAST();
  }
} // namespace transpose