#pragma once
#ifndef TRANSPOSE_H
#define TRANSPOSE_H

#include <cstdint>
#include <stdexcept>
#include "utils/device_context.cuh"
#include "curves/curve_config.cuh"
#include "utils/error_handler.cuh"
#include "utils/utils.h"
#include <sys/time.h>

#define BDIMX 16
#define BDIMY 16
#define IPAD  1

/**
 * @namespace poseidon
 * Implementation of the [Poseidon hash function](https://eprint.iacr.org/2019/458.pdf)
 * Specifically, the optimized [Filecoin version](https://spec.filecoin.io/algorithms/crypto/poseidon/)
 */
namespace transpose {
  inline double seconds()
  {
    struct timeval tp;
    struct timezone tzp;
    int i = gettimeofday(&tp, &tzp);
    return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
  }

  struct TransposeConfig {
    device_context::DeviceContext ctx; /**< Details related to the device such as its id and stream id. */
    bool are_inputs_on_device;  /**< True if inputs are on device and false if they're on host. Default value: false. */
    bool are_outputs_on_device; /**< If true, output is preserved on device, otherwise on host. Default value: false. */
    uint32_t nrows;             /**< Height of the input matrix */
    uint32_t ncols;             /**< Width of the input matrix */
    bool is_async; /**< Whether to run the Poseidon asynchronously. If set to `true`, the poseidon_hash function will be
                    *   non-blocking and you'd need to synchronize it explicitly by running
                    *   `cudaStreamSynchronize` or `cudaDeviceSynchronize`. If set to false, the poseidon_hash
                    *   function will block the current CPU thread. */
  };

  template <typename S>
  class TransposeKernelConfiguration
  {
  public:
    static dim3 get_block() { return dim3(BDIMX, BDIMY); }
    static dim3 get_grid(const TransposeConfig& config)
    {
      int x = (config.ncols + BDIMX - 1) / BDIMX;
      int y = (config.nrows + BDIMY - 1) / BDIMY;
      return dim3((x + 2 - 1) / 2, y);
    }
    static size_t size_bytes(const TransposeConfig& config) { return config.nrows * config.ncols * sizeof(S); }
    static size_t sm_size_bytes() { return (2 * BDIMX + IPAD) * BDIMY * sizeof(S); }
  };
  template <typename S>
  using TKC = TransposeKernelConfiguration<S>;

  template <typename S>
  TransposeConfig default_transpose_config(int nrows, int ncols)
  {
    device_context::DeviceContext ctx = device_context::get_default_device_context();
    TransposeConfig config = {
      ctx,   // ctx
      false, // are_inputes_on_device
      false, // are_outputs_on_device
      nrows, // nrows
      ncols, // ncols
      false  // is_async
    };
    return config;
  }

  /**
   * @param input a pointer to the input data. May be allocated on device or on host, regulated
   * by the config.
   * @param output a pointer to the output data. May be allocated on device or on host, regulated
   * by the config.
   * @param config a struct that encodes various parameters.
   */
  template <typename S>
  cudaError_t transpose_mem(S* input, S* output, const TransposeConfig& config);
} // namespace transpose

#endif