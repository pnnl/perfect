/* -*-Mode: C;-*- */

/**BeginCopyright************************************************************
 *
 * $HeadURL$
 * $Id$
 *
 *---------------------------------------------------------------------------
 * Part of PERFECT Benchmark Suite (hpc.pnnl.gov/projects/PERFECT/)
 *---------------------------------------------------------------------------
 *
 * Copyright ((c)) 2014, Battelle Memorial Institute
 * Copyright ((c)) 2014, Georgia Tech Research Corporation
 * All rights reserved.
 *
 * 1. Battelle Memorial Institute (hereinafter Battelle) and Georgia Tech
 *    Research Corporation (GTRC) hereby grant permission to any person
 *    or entity lawfully obtaining a copy of this software and associated
 *    documentation files (hereinafter "the Software") to redistribute
 *    and use the Software in source and binary forms, with or without
 *    modification.  Such person or entity may use, copy, modify, merge,
 *    publish, distribute, sublicense, and/or sell copies of the
 *    Software, and may permit others to do so, subject to the following
 *    conditions:
 * 
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimers.
 * 
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in
 *      the documentation and/or other materials provided with the
 *      distribution.
 * 
 *    * Other than as used herein, neither the name Battelle Memorial
 *      Institute nor Battelle may be used in any form whatsoever without
 *      the express written consent of Battelle.
 * 
 *      Other than as used herein, neither the name Georgia Tech Research
 *      Corporation nor GTRC may not be used in any form whatsoever
 *      without the express written consent of GTRC.
 * 
 *    * Redistributions of the software in any form, and publications
 *      based on work performed using the software should include the
 *      following citation as a reference:
 * 
 *      Kevin Barker, Thomas Benson, Dan Campbell, David Ediger, Roberto
 *      Gioiosa, Adolfy Hoisie, Darren Kerbyson, Joseph Manzano, Andres
 *      Marquez, Leon Song, Nathan R. Tallent, and Antonino Tumeo.
 *      PERFECT (Power Efficiency Revolution For Embedded Computing
 *      Technologies) Benchmark Suite Manual. Pacific Northwest National
 *      Laboratory and Georgia Tech Research Institute, December 2013.
 *      http://hpc.pnnl.gov/projects/PERFECT/
 *
 * 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
 *    BATTELLE, GTRC, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 *    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **EndCopyright*************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define BLOCK_SIZE 128

__global__ void gpu_bit_reverse_01 (float * w, unsigned int N, unsigned int M, unsigned int bits)
{
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N * M) return;
   //thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(i >= M) return;
   int j_n = blockIdx.y;
   if(j_n >= N) return;
   unsigned int s, shift;
   s = sizeof(int) * CHAR_BIT - 1;
   shift = s - bits + 1;
   unsigned int r;
   float t_real, t_imag;
   r = i;
   unsigned int v = i;
   for (v >>= 1; v; v >>= 1)
   {
      r <<= 1;
      r |= v & 1;
      s--;
   }
   r <<= s;
   r >>= shift;
   if (i < r) {
      t_real = w[(2*j_n*M) + 2 * r];
      t_imag = w[(2*j_n*M) + 2 * r + 1];
      w[(2*j_n*M) + 2 * r] = atomicExch(&w[(2*j_n*M) + 2 * i], t_real);
      w[(2*j_n*M) + 2 * r + 1] = atomicExch(&w[(2*j_n*M) + 2 * i + 1], t_imag);
   }
}

__global__ void gpu_bit_reverse_02 (float * w, unsigned int N, unsigned int M, unsigned int bits)
{

   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * blockDim.y
        + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N * M) return;
   //thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(i >= N) return;
   int i_m = blockIdx.y;
   if(i_m >= M) return;
   unsigned int s, shift;
   s = sizeof(int) * CHAR_BIT - 1;
   shift = s - bits + 1;
   unsigned int r;
   float t_real, t_imag;
   r = i;
   unsigned int v = i;
   for (v >>= 1; v; v >>= 1)
   {
      r <<= 1;
      r |= v & 1;
      s--;
   }
   r <<= s;
   r >>= shift;
   if (i < r) {
      t_real = w[2*(M*r+i_m)    ];
      t_imag = w[2*(M*r+i_m) + 1];
      w[2*(M*r+i_m)    ] = atomicExch(&w[2*(M*i+i_m)    ], t_real);
      w[2*(M*r+i_m) + 1] = atomicExch(&w[2*(M*i+i_m) + 1], t_imag);
   }
}

__global__ void gpu_compute_result_01 (float * data, int N, int M, unsigned int transform_length, float s, float s2)
{
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * blockDim.y
      + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N * M) return;
   int i_n = blockIdx.y;
   if(i_n >= N) return;
   //thread ID for given bock size and number of blocks
   int t_id = (blockIdx.x * blockDim.x) + threadIdx.x;
   int b = t_id * ( 2 * transform_length); 
   if(b >= M) return;
   float w_real;
   float w_imag;
   w_real = 1.0f;
   w_imag = 0.0f;
   for (int a = 0; a < transform_length; a++) {
      int i, j;
      float z_real, z_imag;
      float t_real, t_imag;
      i = b + a;
      j = b + a + transform_length;
      z_real = data[2*i_n*M + 2*j  ];
      z_imag = data[2*i_n*M + 2*j+1];
      t_real = w_real * z_real - w_imag * z_imag;
      t_imag = w_real * z_imag + w_imag * z_real;
      // write the result 
      data[2*i_n*M + 2*j  ]  = data[2*i_n*M + 2*i  ] - t_real;
      data[2*i_n*M + 2*j+1]  = data[2*i_n*M + 2*i+1] - t_imag;
      data[2*i_n*M + 2*i  ] += t_real;
      data[2*i_n*M + 2*i+1] += t_imag;
      t_real = w_real - (s * w_imag + s2 * w_real);
      t_imag = w_imag + (s * w_real - s2 * w_imag);
      w_real = t_real;
      w_imag = t_imag;
   }
}

__global__ void gpu_compute_result_02 (float * data, int N, int M, unsigned int transform_length, float s, float s2)
{
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * blockDim.y
       + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N * M) return;
   int i_m = blockIdx.y;
   if(i_m >= M) return;
   //thread ID for given bock size and number of blocks
   int t_id = (blockIdx.x * blockDim.x) + threadIdx.x;
   int b = t_id * ( 2 * transform_length); 
   if(b >= N) return;
   float w_real;
   float w_imag;
   w_real = 1.0f;
   w_imag = 0.0f;
   for (int a = 0; a < transform_length; a++) {
      int i, j;
      float z_real, z_imag;
      float t_real, t_imag;
      i = b + a;
      j = b + a + transform_length;
      z_real = data[2*(M*j+i_m)    ];
      z_imag = data[2*(M*j+i_m) + 1];
      t_real = w_real * z_real - w_imag * z_imag;
      t_imag = w_real * z_imag + w_imag * z_real;
      // write the result 
      data[2*(M*j+i_m)    ]  = data[2*(M*i+i_m)    ] - t_real;
      data[2*(M*j+i_m) + 1]  = data[2*(M*i+i_m) + 1] - t_imag;
      data[2*(M*i+i_m)    ] += t_real;
      data[2*(M*i+i_m) + 1] += t_imag;
      t_real = w_real - (s * w_imag + s2 * w_real);
      t_imag = w_imag + (s * w_real - s2 * w_imag);
      w_real = t_real;
      w_imag = t_imag;
   }	
}


extern "C" int gpu_fft_2d (float * data, unsigned int N, unsigned int M, unsigned int logn, unsigned int logm, int sign)
{
   //allocate memory for variable data
   float *dev_data;
   CUDA_SAFE(cudaMalloc(&dev_data, 2 * M *N * sizeof(float)));
   CUDA_SAFE(cudaMemcpy(dev_data, data, 2 * M * N * sizeof(float), cudaMemcpyHostToDevice));
   //perform 1D FFT on each row
   int num_blocks = (M + BLOCK_SIZE - 1) / BLOCK_SIZE;
   //compute number of blocks for kernel launch
   dim3 grid_dim_01(num_blocks,N,1);
   dim3 block_dim_01(BLOCK_SIZE);
   gpu_bit_reverse_01<<<grid_dim_01,block_dim_01>>>(dev_data, N, M, logm);
   CUDA_SAFE(cudaDeviceSynchronize());
   // calculation 
   unsigned int transform_length;
   unsigned int bit;
   transform_length = 1;
   for (bit = 0; bit < logn; bit++) {
      float theta, s, t, s2;
      theta = 1.0 * sign * M_PI / (float) transform_length;
      s = sin (theta);
      t = sin (0.5 * theta);
      s2 = 2.0 * t * t;
      int num_blocks = (M / (2 * transform_length) + BLOCK_SIZE - 1) / BLOCK_SIZE;
      //compute number of blocks for kernel launch
      dim3 grid_dim_01(num_blocks,N,1);
      dim3 block_dim_01(BLOCK_SIZE);
      gpu_compute_result_01<<<grid_dim_01,block_dim_01>>>(dev_data, N, M, transform_length, s, s2);
      CUDA_SAFE(cudaDeviceSynchronize());
      transform_length *= 2;
   }
   //end of 1D FFT on each row
   //perform 1D FFT on each column
   num_blocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
   //compute number of blocks for kernel launch
   dim3 grid_dim_02(num_blocks,M,1);
   dim3 block_dim_02(BLOCK_SIZE);
   gpu_bit_reverse_02<<<grid_dim_02,block_dim_02>>>(dev_data, N, M, logm);
   CUDA_SAFE(cudaDeviceSynchronize());
   transform_length = 1;
   for (bit = 0; bit < logn; bit++) {
      float theta, s, t, s2;
      theta = 1.0 * sign * M_PI / (float) transform_length;
      s = sin (theta);
      t = sin (0.5 * theta);
      s2 = 2.0 * t * t;
      int num_blocks = (N / (2 * transform_length) + BLOCK_SIZE - 1) / BLOCK_SIZE;
      //compute number of blocks for kernel launch
      dim3 grid_dim_02(num_blocks,M,1);
      dim3 block_dim_02(BLOCK_SIZE);
      gpu_compute_result_02<<<grid_dim_02,block_dim_02>>>(dev_data, N, M, transform_length, s, s2);
      CUDA_SAFE(cudaDeviceSynchronize());
      transform_length *= 2;
   }
   //end of 1D FFT on each column
   //copy variable data back from Device to Host
   CUDA_SAFE(cudaMemcpy(data, dev_data, (2 * M *N * sizeof(float)), cudaMemcpyDeviceToHost));
   //free device variables
   CUDA_SAFE(cudaFree(dev_data));
   return 0;
}

