#define CURVE_ID 5
#define NROWS    4096
#define NCOLS    4096
#define DATA_ON_DEVICE
// #define DEBUG
#include "../../curves/curve_config.cuh"
#include "../../utils/device_context.cuh"
#include "transpose.cu"

#ifndef __CUDA_ARCH__
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>

using namespace transpose;
using namespace curve_config;

#define START_TIMER(timer) auto timer##_start = std::chrono::high_resolution_clock::now();
#define END_TIMER(timer, msg)                                                                                          \
  printf("%s: %.0f ms\n", msg, FpMilliseconds(std::chrono::high_resolution_clock::now() - timer##_start).count());
#define END_TIMER_BND(timer, size, msg)                                                                                \
  {                                                                                                                    \
    auto timer##_end = std::chrono::high_resolution_clock::now();                                                      \
    float timer##_duration = FpMilliseconds(timer##_end - timer##_start).count();                                      \
    printf(                                                                                                            \
      "%s: %f ms, %f GB\n", msg, timer##_duration, (float)size / 1024 / 1024 / 1024 * 2 / timer##_duration * 1000);  \
  }

int main(int argc, char* argv[])
{
  using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
  using FpMicroseconds = std::chrono::duration<float, std::chrono::microseconds::period>;

  START_TIMER(allocation_timer);
  // Prepare input data of [0, 1, 2 ... (number_of_blocks * arity) - 1]
  int nrows = argc > 1 ? 1 << atoi(argv[1]) : NROWS;
  int ncols = argc > 2 ? 1 << atoi(argv[2]) : NCOLS;
  printf("nrows = %d, ncols = %d\n", nrows, ncols);
  uint32_t size_bytes = nrows * ncols * sizeof(scalar_t);
  scalar_t input = scalar_t::zero();
  scalar_t* idata = static_cast<scalar_t*>(malloc(size_bytes));
  for (uint32_t i = 0; i < nrows * ncols; i++) {
    idata[i] = input;
    input = input + scalar_t::one();
  }
#ifdef DATA_ON_DEVICE
  scalar_t* in_ptr;
  CHK_IF_RETURN(cudaMalloc(&in_ptr, size_bytes))
  CHK_IF_RETURN(cudaMemcpy(in_ptr, idata, size_bytes, cudaMemcpyHostToDevice))
#endif
  END_TIMER(allocation_timer, "Allocate mem and fill input");

  scalar_t* odata = static_cast<scalar_t*>(malloc(size_bytes));
#ifdef DATA_ON_DEVICE
  scalar_t* out_ptr;
  CHK_IF_RETURN(cudaMalloc(&out_ptr, size_bytes))
  CHK_IF_RETURN(cudaMemset(out_ptr, 0, size_bytes))
#endif

  printf("TransposeConfig\n");
  TransposeConfig config = default_transpose_config<scalar_t>(nrows, ncols);
#ifdef DATA_ON_DEVICE
  config.are_inputs_on_device = true;
  config.are_outputs_on_device = true;
#endif
  printf("TransposeConfig done\n");

// warm up
#ifdef DATA_ON_DEVICE
  transpose_mem<curve_config::scalar_t>(in_ptr, out_ptr, config);
#else
  transpose_mem<curve_config::scalar_t>(idata, odata, config);
#endif

  CHK_IF_RETURN(cudaDeviceSynchronize())
  printf("transpose_mem\n");
  START_TIMER(transpose_timer);
#ifdef DATA_ON_DEVICE
  transpose_mem<curve_config::scalar_t>(in_ptr, out_ptr, config);
#else
  transpose_mem<curve_config::scalar_t>(idata, odata, config);
#endif
  END_TIMER_BND(transpose_timer, size_bytes, "Transpose")
  printf("transpose_mem done\n");
#ifdef DATA_ON_DEVICE
  CHK_IF_RETURN(cudaMemcpy(odata, out_ptr, size_bytes, cudaMemcpyDeviceToHost))
#endif

#ifdef DEBUG
  for (int i = 0; i < nrows; i++) {
    for (int j = 0; j < ncols; j++) {
      std::cout << idata[i * ncols + j] << ", ";
    }
    printf("\n");
  }
  printf("==============================================================\n");

  for (int i = 0; i < ncols; i++) {
    for (int j = 0; j < nrows; j++) {
      std::cout << odata[i * nrows + j] << ", ";
    }
    printf("\n");
  }
#endif

  for (int i = 0; i < ncols; i++) {
    for (int j = 0; j < nrows; j++) {
      assert(odata[i * nrows + j] == idata[j * ncols + i]);
    }
  }

#ifdef DATA_ON_DEVICE
  CHK_IF_RETURN(cudaFree(in_ptr));
  CHK_IF_RETURN(cudaFree(out_ptr));
#endif
  free(idata);
  free(odata);
}

#endif