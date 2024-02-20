// #define DEBUG

#define CURVE_ID 5
#include "../../curves/curve_config.cuh"
#include "../../utils/device_context.cuh"
#include "poseidon.cu"

#ifndef __CUDA_ARCH__
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>

using namespace poseidon;
using namespace curve_config;

#define A 11
#define T (A + 1)

#define START_TIMER(timer) auto timer##_start = std::chrono::high_resolution_clock::now();
#define END_TIMER(timer, msg)                                                                                          \
  printf("%s: %.0f ms\n", msg, FpMilliseconds(std::chrono::high_resolution_clock::now() - timer##_start).count());

int main(int argc, char* argv[])
{
  using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
  using FpMicroseconds = std::chrono::duration<float, std::chrono::microseconds::period>;

  // Load poseidon constants
  START_TIMER(timer_const);
  device_context::DeviceContext ctx = device_context::get_default_device_context();
  PoseidonConstants<scalar_t> constants;
  printf("init_optimized_poseidon_constants\n");
  init_optimized_poseidon_constants<scalar_t>(A, ctx, &constants);
  printf("init_optimized_poseidon_constants done\n");
  END_TIMER(timer_const, "Load poseidon constants");

  START_TIMER(allocation_timer);
  // Prepare input data of [0, 1, 2 ... (number_of_blocks * arity) - 1]
  int number_of_blocks = argc > 1 ? 1 << atoi(argv[1]) : 1;
  printf("number_of_blocks = %d\n", number_of_blocks);
  scalar_t input = scalar_t::zero();
  scalar_t* in_ptr = static_cast<scalar_t*>(malloc(number_of_blocks * T * sizeof(scalar_t)));
  for (uint32_t i = 0; i < number_of_blocks * T; i++) {
    in_ptr[i] = input;
    input = input + scalar_t::one();
  }
  END_TIMER(allocation_timer, "Allocate mem and fill input");

  printf("out_ptr\n");
  scalar_t* out_ptr = static_cast<scalar_t*>(malloc(number_of_blocks * 4 * sizeof(scalar_t)));
  printf("out_ptr done\n");

  START_TIMER(poseidon_timer);
  printf("PoseidonConfig\n");
  PoseidonConfig config = default_poseidon_config<scalar_t>(T);
  printf("PoseidonConfig done\n");
  printf("poseidon_hash\n");
  poseidon_hash<curve_config::scalar_t, T>(in_ptr, out_ptr, number_of_blocks, constants, config);
  printf("poseidon_hash done\n");
  END_TIMER(poseidon_timer, "Poseidon")

  scalar_t expected[4] = {
    {0xfc5b8e9e, 0xd64e1e3e}, {0x020aaa47, 0x53666633}, {0x7c6a8825, 0xd4028559}, {0xe81231d2, 0x613a4f81}};

  if (number_of_blocks == 1) {
    for (int i = 0; i < 4; i++) {
#ifdef DEBUG
      std::cout << out_ptr[i] << std::endl;
#endif
      assert((out_ptr[i] == expected[i]));
    }
    printf("Expected output matches\n");
  }

  free(in_ptr);
  free(out_ptr);
}

#endif