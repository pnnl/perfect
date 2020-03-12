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

#define BLOCK_SIZE_A 32
#define SM_ARRAY 512
#define BLOCK_SIZE_B 64
#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

__global__ void gpu_n_squared_sort (float *value_, int *index_, int len)
{
   __shared__ float value[SM_ARRAY];
   __shared__ int index[SM_ARRAY];
   int num_blocks = ( len + SM_ARRAY - 1) / SM_ARRAY;
   if (blockIdx.x >= num_blocks) return;
   //copy data to shared memory
   for (int i = 0; i < (SM_ARRAY/BLOCK_SIZE_A); i++)
   {
      value[threadIdx.x] = value_[blockIdx.x * SM_ARRAY + i*BLOCK_SIZE_A + threadIdx.x];
      index[threadIdx.x] = index_[blockIdx.x * SM_ARRAY + i*BLOCK_SIZE_A + threadIdx.x];
   }
   for (int jj = 0; jj < (SM_ARRAY + BLOCK_SIZE_A - 1)/BLOCK_SIZE_A; jj++)
   {
      int j = jj + threadIdx.x * (int)((SM_ARRAY + BLOCK_SIZE_A - 1)/BLOCK_SIZE_A);
      if (value[j] > value[j+1] && j < (SM_ARRAY-1) )
      {
         double val_tmp;
         int idx_tmp;
         val_tmp = value[j];
         value[j] = value[j+1];
         value[j+1] = val_tmp;
         idx_tmp = index[j];
         index[j] = index[j+1];
         index[j+1] = idx_tmp;
      }
   }
   __syncthreads();
   //copy data back from shared memory
   for (int i = 0; i < (SM_ARRAY/BLOCK_SIZE_A); i++)
   {
      value_[blockIdx.x * SM_ARRAY + i*BLOCK_SIZE_A + threadIdx.x] = value[threadIdx.x];
      index_[blockIdx.x * SM_ARRAY + i*BLOCK_SIZE_A + threadIdx.x] = index[threadIdx.x];
   }
}

extern "C" void gpu_sort(float * value, int * index, int len, int gpu_index)
{
   float *dev_value;
   int *dev_index;

   CUDA_SAFE_MALLOC(dev_value, sizeof(float) * len);
   CUDA_SAFE_MALLOC(dev_index, sizeof(int) * len);

   CUDA_SAFE(cudaMemcpy(dev_value, value, sizeof(float) * len, cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_index, index, sizeof(int) * len, cudaMemcpyHostToDevice));

   if(gpu_index == 1)
   {
      //compute number of blocks for kernel launch
      int num_blocks = ( len + SM_ARRAY - 1) / SM_ARRAY;

      //compute number of blocks for kernel launch
      dim3 grid_dim_a(num_blocks,(len-1),1);
      dim3 block_dim_a(BLOCK_SIZE_A);

      gpu_n_squared_sort<<<grid_dim_a,block_dim_a,(len/BLOCK_SIZE_A)>>>(dev_value, dev_index, len);
      //synchronize threads
      CUDA_SAFE(cudaDeviceSynchronize());
   }
   else if(gpu_index == 2)
   {
      /* currently not implemented as it is too slow*/
   }
   //copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(value, dev_value, sizeof(float) * len, cudaMemcpyDeviceToHost));
   //free device variables
   cudaFree(dev_value);
   cudaFree(dev_index);
}
