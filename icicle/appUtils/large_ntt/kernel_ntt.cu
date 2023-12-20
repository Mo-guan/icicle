#ifndef K_NTT
#define K_NTT
#pragma once

/*

Copyright (c) 2023 Yrrid Software, Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the �Software�), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions
of the Software.

THE SOFTWARE IS PROVIDED �AS IS�, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

*/

// #include <stdio.h>
// #include <stdint.h>
// #include "asm.cu"
// #include "Sampled96.cu"

// typedef Sampled96 Math;

#include "thread_ntt.cu"

// EXCHANGE_OFFSET = 129*64*8

#define DATA_OFFSET 0

__launch_bounds__(384)
__global__ void ntt1024(uint64_t* out, uint64_t* in, uint32_t* next, uint32_t count) {
  NTTEngine32 engine;
  uint32_t    dataIndex=0;

  #ifdef COMPUTE_ONLY
    bool        first=true;
  #endif
  
  engine.initializeRoot();
    
  while(true) {
    if((threadIdx.x & 0x1F)==0)
      dataIndex=atomicAdd(next, 1);
    dataIndex=__shfl_sync(0xFFFFFFFF, dataIndex, 0);  //send value to all threads in warp    
    if(dataIndex<count) {
      #if defined(COMPUTE_ONLY)
        if(first)
          engine.loadGlobalData(in, dataIndex);
        first=false;
      #else
        engine.loadGlobalData(in, dataIndex);
      #endif
    }
    else {
      if(dataIndex==count + (gridDim.x*blockDim.x>>5) - 1) { //didn't understand condition
        // last one to finish, reset the counter
        atomicExch(next, 0);
      }
      return;
    }
    #pragma unroll 1
    for(uint32_t phase=0;phase<2;phase++) {
      // ntt32 produces a lot of instructions, so we put this in a loop
      engine.ntt32(); 
      if(phase==0) {
        engine.storeSharedData(DATA_OFFSET);
        __syncwarp(0xFFFFFFFF);
        engine.loadSharedDataAndTwiddle32x32(DATA_OFFSET);
      }
    }
    engine.storeGlobalData(out, dataIndex);
  }
}


__global__ void thread_ntt_kernel(test_scalar* out, test_scalar* in, uint32_t* next, uint32_t count) {
  NTTEngine engine;
  uint32_t    dataIndex=blockIdx.x*blockDim.x+threadIdx.x;

  #ifdef COMPUTE_ONLY
    bool        first=true;
  #endif
  
  // engine.initializeRoot();
  engine.loadGlobalDataDep(in, dataIndex);
  // engine.ntt4_4();
  // for (int i = 0; i < 16; i++)
  // {
  //   engine.X[i] = in[i];
  // }
  

  // for (int i = 0; i < 100; i++)
  // {
    // engine.ntt16win_lowreg();
    // engine.ntt8_2();
    // engine.ntt8_2();
    // engine.ntt16();
    // engine.ntt16_win8ct2();
    engine.ntt16win();
    // engine.X[2] = engine.X[2]*engine.X[0];
  // }
  // engine.ntt16();
  // engine.X[0] = engine.X[0] + engine.X[1];
  // out[0] = out[0] + out[1];
  // out[0] = test_scalar::zero();
  // engine.storeGlobalData(out, dataIndex);
  // engine.storeGlobalData8_2(out, dataIndex);
  // engine.storeGlobalData16(out, dataIndex);
    
  // while(true) {
  //   if((threadIdx.x & 0x1F)==0)
  //     dataIndex=atomicAdd(next, 1);
  //   dataIndex=__shfl_sync(0xFFFFFFFF, dataIndex, 0);      
  //   if(dataIndex<count) {
  //     #if defined(COMPUTE_ONLY)
  //       if(first)
  //         engine.loadGlobalData(in, dataIndex);
  //       first=false;
  //     #else
  //       engine.loadGlobalData(in, dataIndex);
  //     #endif
  //   }
  //   else {
  //     if(dataIndex==count + (gridDim.x*blockDim.x>>5) - 1) {
  //       // last one to finish, reset the counter
  //       atomicExch(next, 0);
  //     }
  //     return;
  //   }
  //   #pragma unroll 1
  //   for(uint32_t phase=0;phase<2;phase++) {
  //     // ntt32 produces a lot of instructions, so we put this in a loop
  //     engine.ntt32(); 
  //     if(phase==0) {
  //       engine.storeSharedData(DATA_OFFSET);
  //       __syncwarp(0xFFFFFFFF);
  //       engine.loadSharedDataAndTwiddle32x32(DATA_OFFSET);
  //     }
  //   }
    engine.storeGlobalDataDep(out, dataIndex);
  // }
}

__launch_bounds__(64)
// __global__ void ntt_kernel_split_transpose(test_scalar* out, test_scalar* in) {
__global__ void ntt_kernel_split_transpose(uint4* out, uint4* in) {
  NTTEngine engine;
  uint32_t    dataIndex=blockIdx.x*blockDim.x+threadIdx.x;

  // if (blockIdx.x !=1) return;

  // if (blockIdx.x ==0 && threadIdx.x ==0) printf("start kernel\n");  
  // __shared__ uint4 shmem[2048*3];
  extern __shared__ uint4 shmem[];
  // if (blockIdx.x ==0 && threadIdx.x ==0) printf("shmem\n");

  // #ifdef COMPUTE_ONLY
  //   bool        first=true;
  // #endif
  
  // engine.initializeRoot();
  // engine.loadGlobalData(in, dataIndex);
  // engine.ntt4_4();
  // for (int i = 0; i < 100000; i++)
  // {
    // engine.ntt16win();
    // engine.ntt16win_lowreg();
    // engine.ntt8_2();
    // engine.ntt8_2();
    // engine.ntt16();
  // engine.ntt16_win8ct2();
  // }
  // engine.ntt16();
  // engine.X[0] = engine.X[0] + engine.X[1];
  // out[0] = out[0] + out[1];
  // out[0] = test_scalar::zero();
  // engine.storeGlobalData(out, dataIndex);
  // engine.storeGlobalData8_2(out, dataIndex);
  // engine.storeGlobalData16(out, dataIndex);
    
  // while(true) {
    // if((threadIdx.x & 0x1F)==0)
    //   dataIndex=atomicAdd(next, 1);
    // dataIndex=__shfl_sync(0xFFFFFFFF, dataIndex, 0);      
    // if(dataIndex<count) {
    //   #if defined(COMPUTE_ONLY)
    //     if(first)
    //       engine.loadGlobalData(in, dataIndex);
    //     first=false;
    //   #else
    //     engine.loadGlobalData(in, dataIndex);
    //   #endif
    // }
    // else {
    //   if(dataIndex==count + (gridDim.x*blockDim.x>>5) - 1) {
    //     // last one to finish, reset the counter
    //     atomicExch(next, 0);
    //   }
    //   return;
    // }
    // if (threadIdx.x!=0) return;
    // engine.loadGlobalDataDep(in, dataIndex);
    // engine.loadGlobalData(in,blockIdx.x*512*2,1,64*8); //todo - change function to fit global ntt
    // engine.loadGlobalData(in,blockIdx.x*512*2,1,256*8); //todo - change function to fit global ntt
    // __syncthreads();
    // if (blockIdx.x ==0 && threadIdx.x ==0) printf("load global\n");
    // engine.externalTwiddles(); //todo
    // engine.twiddles256();
    // engine.ntt16_win8ct2();
    // engine.twiddles256();
    // engine.ntt16_win8ct2();
    // #pragma unroll 1
    // for (uint32_t i=0;i<100;i++) {
    #pragma unroll 1
    for (uint32_t phase=0;phase<2;phase++) {
      // ntt32 produces a lot of instructions, so we put this in a loop
      // engine.ntt16_win8ct2();
      // engine.plus();
      // engine.twiddles256();
      // engine.load_twiddles(threadIdx.x&0x7);
      // engine.twiddles64(threadIdx.x&0x7);
      // engine.ntt16_win8ct2();
      // engine.ntt8win();
      // engine.ntt16_win8ct2();
      // engine.twiddles256();
      // engine.ntt16_win8ct2();
      // engine.ntt16win();
      // engine.twiddles256();
      // engine.ntt16();
      // if(phase==0) {
      //   engine.SharedDataColumns2(shmem, true, false); //store low
      //   __syncthreads();
      //   // if (blockIdx.x ==0 && threadIdx.x ==0) printf("store shmem low\n");
      //   // if (blockIdx.x ==0 && threadIdx.x ==0){
      //   //   for (int i = 0; i < 512; i++)
      //   //   {
      //   //     if (i%32==0) printf("\n");
      //   //     if (i%256==0) printf("\n");
      //   //     printf("%d, ",shmem[i].w);
      //   //   }
      //   // }
      //   // __syncthreads();
      //   engine.SharedDataRows2(shmem, false, false); //load low
      //   // if (blockIdx.x ==0 && threadIdx.x ==0) printf("load shmem low\n");
      //   // __syncthreads(); //can avoid with switching rows and columns
      //   // if (blockIdx.x ==0 && threadIdx.x ==1){
      //   //   for (int i = 0; i < 16; i++)
      //   //   {
      //   //     printf("\n");
      //   //     printf("%d, ",engine.X[i].limbs_storage.limbs[0]);
      //   //   }
      //   // }
      //   engine.SharedDataRows2(shmem, true, true); //store high
      //   __syncthreads();
      //   // if (blockIdx.x ==0 && threadIdx.x ==0) printf("store shmem high\n");
      //   // if (blockIdx.x ==0 && threadIdx.x ==0){
      //   //   for (int i = 0; i < 2048; i++)
      //   //   {
      //   //     if (i%16==0) printf("\n");
      //   //     if (i%256==0) printf("\n");
      //   //     printf("%d, ",shmem[i].w);
      //   //   }
      //   // }
      //   // __syncthreads();
      //   engine.SharedDataColumns2(shmem, false, true); //load high
      //   // if (blockIdx.x ==0 && threadIdx.x ==0) printf("load shmem high\n");
      //   // engine.twiddles256();
      //   // engine.ntt16_win8ct2();
      //   // engine.twiddles256();
      // }
    // }
    }
    // #pragma unroll 1
    // for (uint32_t i=0;i<100*2;i++) {
    // // #pragma unroll 1
    // // for (uint32_t phase=0;phase<2;phase++) {
    //   // ntt32 produces a lot of instructions, so we put this in a loop
    //   // engine.ntt16_win8ct2();
    //   // if (i%2) engine.twiddles256();
    //   engine.twiddles256();
    //   engine.ntt16_win8ct2();
    //   // engine.ntt16win();
    // }
    // engine.storeGlobalData(out,blockIdx.x*512*2,1,64*8); //todo - change function to fit global ntt
    // engine.storeGlobalData(out,blockIdx.x*512*2,1,256*8); //todo - change function to fit global ntt
    // engine.storeGlobalDataDep(out, dataIndex); //todo - change function to fit global ntt
  // }
}

__launch_bounds__(64)
__global__ void ntt64(uint4* in, uint4* out, uint4* twiddles, uint4* internal_twiddles, uint32_t log_size, uint32_t tw_log_size, uint32_t data_stride, uint32_t twiddle_stride) {
// __global__ void ntt64(uint4* in, uint4* out, uint4* twiddles, uint32_t log_size, uint32_t tw_log_size, uint32_t data_stride, uint32_t twiddle_stride) {
  NTTEngine engine;
  extern __shared__ uint4 shmem[];

  // if (twiddle_stride && threadIdx.x>15) return;
  
  // engine.initializeRoot(data_stride>1);
  engine.loadInternalTwiddles(internal_twiddles, data_stride>1);
    
  // #pragma unroll 1
  // for (int i = 0; i < 260; i++) //todo - function of size
  // {
  engine.loadGlobalData(in, data_stride, twiddle_stride, log_size);
  if (twiddle_stride) {
    engine.loadExternalTwiddles(twiddles, data_stride, twiddle_stride, log_size, tw_log_size);
    engine.twiddlesExternal();
  }


    // if (twiddle_stride){
    //   engine.loadGlobalData(in,blockIdx.x*64*8,data_stride,1<<log_size); //todo - parametize
    //   engine.loadExternalTwiddles(twiddles,blockIdx.x*64*8,twiddle_stride,1<<tw_log_size); //todo - parametize
    //   engine.twiddlesExternal();
    // }
    // else{
    //   engine.loadGlobalData(in,blockIdx.x*8,data_stride,1<<log_size); //todo - parametize
    // }
    // engine.twiddles64();

    #pragma unroll 1
    for(uint32_t phase=0;phase<2;phase++) {
      // this code produces a lot of instructions, so we put this in a loop
      engine.ntt8win(); 
      if(phase==0) {
        engine.SharedDataColumns2(shmem, true, false, data_stride>1); //store low
        __syncthreads();
        engine.SharedDataRows2(shmem, false, false, data_stride>1); //load low
        engine.SharedDataRows2(shmem, true, true, data_stride>1); //store high
        __syncthreads();
        engine.SharedDataColumns2(shmem, false, true, data_stride>1); //load high
        engine.twiddles64();
      }
    }

    engine.storeGlobalData(twiddle_stride? out :in, data_stride, log_size);
  // }
}

// __global__ void generate_base_tables(test_scalar* w6_table, test_scalar* w12_table, test_scalar* w18_table){

// test_scalar w, t;

// switch (threadIdx.x)
// {
// case 0:
//   w = test_scalar::omega(6);
//   t = test_scalar::one();
//   for (int i = 0; i < 64; i++)
//   {
//     // w6_table[i] = t.load_half(false);
//     // w6_table[i + 64] = t.load_half(true);
//     w6_table[i] = t;
//     t = t * w;
//   }
//   break;
// case 1:
//   w = test_scalar::omega(12);
//   t = test_scalar::one();
//   for (int i = 0; i < 64; i++)
//   {
//     w12_table[i] = t;
//     t = t * w;
//   }
//   break;
// case 2:
//   w = test_scalar::omega(18);
//   t = test_scalar::one();
//   for (int i = 0; i < 64; i++)
//   {
//     w18_table[i] = t;
//     t = t * w;
//   }
//   break;

// default:
//   break;
// }

// }

__global__ void generate_base_table(int x, uint4* base_table){
  
  test_scalar w = test_scalar::omega(x);
  test_scalar t = test_scalar::one();
  for (int i = 0; i < 64; i++)
  {
    // base_table[i] = t;
    base_table[i] = t.load_half(false);
    base_table[i+64] = t.load_half(true);
    t = t * w;
  }
}

__global__ void generate_twiddle_combinations(uint4* w6_table, uint4* w12_table, uint4* w18_table, uint4* w24_table, uint4* twiddles){
  
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  int exp = (tid & 0x3f) * (tid >> 6);
  test_scalar w6,w12,w18,w24;
  w6.store_half(w6_table[exp >> 18], false);
  w6.store_half(w6_table[(exp >> 18) + 64], true);
  w12.store_half(w12_table[((exp >> 12) & 0x3f)], false);
  w12.store_half(w12_table[((exp >> 12) & 0x3f) + 64], true);
  w18.store_half(w18_table[((exp >> 6) & 0x3f)], false);
  w18.store_half(w18_table[((exp >> 6) & 0x3f) + 64], true);
  w24.store_half(w24_table[(exp & 0x3f)], false);
  w24.store_half(w24_table[(exp & 0x3f) + 64], true);
  // test_scalar t = w6_table[exp >> 18] * w12_table[(exp >> 12) & 0x3f] * w18_table[(exp >> 6) & 0x3f] * w24_table[exp & 0x3f];
  test_scalar t = w6 * w12 * w18 * w24;
  // test_scalar t = w6_table[exp >> 12] * w12_table[(exp >> 6) & 0x3f] * w18_table[exp & 0x3f];
  // test_scalar t = w6_table[exp >> 6] * w12_table[exp & 0x3f];
  twiddles[tid] = t.load_half(false);
  twiddles[tid + (1<<24)] = t.load_half(true);
  // twiddles[tid + (1<<18)] = t.load_half(true);
  // twiddles[tid + (1<<12)] = t.load_half(true);

}
// __global__ void generate_internal_twiddles(uint4* twiddles){
//   test_scalar w = test_scalar::omega(6);
//   test_scalar t = test_scalar::one();
//   for (int i = 0; i < 64; i++)
//   {
//     twiddles[i] = t.load_half(false);
//     twiddles[i+64] = t.load_half(true);
//     t = t * w;
//   }
// }

uint4* generate_external_twiddles(uint4* twiddles, uint32_t log_size){

  uint4* w6_table;
  uint4* w12_table;
  uint4* w18_table;
  uint4* w24_table;
  cudaMalloc((void**)&w6_table, sizeof(uint4)*64*2);
  cudaMalloc((void**)&w12_table, sizeof(uint4)*64*2);
  cudaMalloc((void**)&w18_table, sizeof(uint4)*64*2);
  cudaMalloc((void**)&w24_table, sizeof(uint4)*64*2);
  generate_base_table<<<1,1>>>(6, w6_table);
  generate_base_table<<<1,1>>>(12, w12_table);
  generate_base_table<<<1,1>>>(18, w18_table);
  generate_base_table<<<1,1>>>(24, w24_table);
  generate_twiddle_combinations<<<1024*64,256>>>(w6_table, w12_table, w18_table, w24_table, twiddles);
  // generate_twiddle_combinations<<<1024,256>>>(w6_table, w12_table, w18_table, twiddles);
  // generate_twiddle_combinations<<<16,256>>>(w6_table, w12_table, w18_table, twiddles);
  return w6_table;
}

void new_ntt(uint4* in, uint4* out, uint4* twiddles, uint4* internal_twiddles, uint32_t log_size, uint32_t tw_log_size){
// void new_ntt(uint4* in, uint4* out, uint4* twiddles, uint32_t log_size){
  uint32_t nof_stages = 2;
  // for (int i = 0; i < nof_stages; i++)
  // {
  //   ntt64<<<1<<(log_size-9),64,8*64*sizeof(uint4)>>>(in, out, twiddles, log_size, i? 1 : 64, i? 1 : 0);
  // }
  // ntt64<<<8,64,8*64*sizeof(uint4)>>>(in, out, twiddles, internal_twiddles, log_size,24, 64,0);
  // ntt64<<<8,64,8*64*sizeof(uint4)>>>(in, out, twiddles, internal_twiddles, log_size,24, 1, 1);
  ntt64<<<8*64,64,8*64*sizeof(uint4)>>>(in, out, twiddles, internal_twiddles, log_size, tw_log_size, 64*64,0);
  ntt64<<<8*64,64,8*64*sizeof(uint4)>>>(in, out, twiddles, internal_twiddles, log_size, tw_log_size, 64, 64*64);
  ntt64<<<8*64,64,8*64*sizeof(uint4)>>>(in, out, twiddles, internal_twiddles, log_size, tw_log_size, 1, 1);
}

#endif