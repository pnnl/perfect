/**************************/
/***    UNCLASSIFIED    ***/
/**************************/

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


/***

  ALL SOURCE CODE PRESENT IN THIS FILE IS UNCLASSIFIED AND IS
  BEING PROVIDED IN SUPPORT OF THE DARPA PERFECT PROGRAM.

  THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY, EXPRESSED, IMPLIED, 
  OR OTHERWISE INFERRED. USE AND SUITABILITY FOR ANY PARTICULAR
  APPLICATION IS SOLELY THE RESPONSIBILITY OF THE IMPLEMENTER. 
  NO CLAIM OF SUITABILITY FOR ANY APPLICATION IS MADE.
  USE OF THIS CODE FOR ANY APPLICATION RELEASES THE AUTHOR
  AND THE US GOVT OF ANY AND ALL LIABILITY.

  THE POINT OF CONTACT FOR QUESTIONS REGARDING THIS SOFTWARE IS:

  US ARMY RDECOM CERDEC NVESD, RDER-NVS-SI (JOHN HODAPP), 
  10221 BURBECK RD, FORT BELVOIR, VA 22060-5806

  THIS HEADER SHALL REMAIN PART OF ALL SOURCE CODE FILES.

 ***/

/*
   This CUDA code is implementation of the kernel - histogram equalization
   as per C code in file histeq.c provided with header as above

   The cuda implementation done by the team at PNNL
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../inc/histeq.h"

/*
 define block size for kernel launch, could be variable
 tested from a value of 64 upto a maximum of 256
 use of shared memory limit the maximum block size
 */

#define BLOCK_SIZE 256
#define MAX_NFILTERROWS 10
#define MAX_NFILTERCOLS 10

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) {    \
   printf("CUDA CALL FAILED AT %d\n", __LINE__ );   \
   exit(1);}

#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

//define device variable
__device__ int CDFblock = 0;

__global__ void gpu_hist(int *streamA, int *h, int nRows, int nCols, int nBins)
{
   //shared variables using shared memory
   __shared__ int shared_streamA[BLOCK_SIZE];
   __shared__ int shared_counter[BLOCK_SIZE];
   //thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= nRows * nCols ) return;
   //block level reduction used to computer array variable h
   shared_streamA[threadIdx.x] = streamA[i];
   shared_counter[threadIdx.x] = 1;
   __syncthreads();
   for(int offset = blockDim.x>>1; offset > 0; offset >>=1)
   {
      if(threadIdx.x < offset)
      {
         if(shared_streamA[threadIdx.x] == 
            shared_streamA[threadIdx.x + offset])
         {
            shared_counter[threadIdx.x] += shared_counter[threadIdx.x + offset];
            shared_counter[threadIdx.x + offset] = 0;
         }
      }
      __syncthreads();
   }
   if(shared_counter[threadIdx.x] > 0)
   {
      if (shared_streamA[threadIdx.x] >= nBins)
      {
         printf("Line %d, Range Error in hist() -- using max val ---", 
                __LINE__);
         atomicAdd(&h[nBins-1], shared_counter[threadIdx.x]);
      }
      else
      {
         atomicAdd(&h[(int)shared_streamA[threadIdx.x]], 
                   shared_counter[threadIdx.x]);
      }
   }
}

__global__ void gpu_histEq_CDF_01(int *h, double *CDF, double *CDF_BLOCK, 
                                  int nInpBins)
{
   /* reduction method used to compute variable CDF
      started in kernel gpu_histEq_CDF_01 and
      completed in kernel gpu_histEq_CDF_02
   */
   //shared variables using shared memory
   __shared__ int shared_h[BLOCK_SIZE];
   __shared__ double shared_CDF[BLOCK_SIZE];
   //thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= nInpBins ) return;
   int j;
   shared_h[threadIdx.x] = h[i];
   shared_CDF[threadIdx.x] = 0.e0;
   __syncthreads();
   for (j = 0; j < BLOCK_SIZE; j++)
   {
      if(threadIdx.x >= j) shared_CDF[threadIdx.x] += (double) shared_h[j];
   }
   __syncthreads();
   CDF[i] = shared_CDF[threadIdx.x];
   if(threadIdx.x == 0) 
   {
      CDF_BLOCK[blockIdx.x] = shared_CDF[blockDim.x - 1];
   }
}

__global__ void gpu_histEq_CDF_02(double *CDF, double *CDF_BLOCK, int nInpBins)
{
/* reduction method used to compute variable CDF
   started in kernel gpu_histEq_CDF_01 and
   completed in kernel gpu_histEq_CDF_02
*/
   //shared variables using shared memory
   __shared__ double shared_CDF_BLOCK[BLOCK_SIZE];
   __shared__ double shared_CDF[BLOCK_SIZE];
   //thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i < blockDim.x ) return;
   if ( i >= nInpBins ) return;
   shared_CDF[threadIdx.x] = CDF[i];
   shared_CDF_BLOCK[threadIdx.x] = CDF[i];
   if(threadIdx.x == 0)
   { 
      for (int j = 0; j < blockIdx.x; j++ )
      {
         shared_CDF_BLOCK[0] += CDF_BLOCK[j];
      }
   }
   __syncthreads();
   shared_CDF[threadIdx.x] += shared_CDF_BLOCK[0];
   CDF[i] = shared_CDF[threadIdx.x];
}

__global__ void gpu_histEq_CDFmin_01(double *CDF, double * CDF_BLOCK, 
                                     int nInpBins)
{
/* reduction method used to compute variable CDFmin
   started in kernel gpu_histEq_CDFmin_01 and
   completed in kernel gpu_histEq_CDFmin_02
*/
   //shared variables using shared memory
   __shared__ double shared_CDF[BLOCK_SIZE];
   // thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= nInpBins ) return;
   shared_CDF[threadIdx.x] = CDF[i];
   __syncthreads();
   for(int offset = blockDim.x>>1; offset > 0; offset >>=1)
   {
      if(threadIdx.x < offset)
      {
         if(shared_CDF[threadIdx.x] > shared_CDF[threadIdx.x + offset])
         {
            shared_CDF[threadIdx.x] = shared_CDF[threadIdx.x+offset];
         }
      }
      __syncthreads();
   }
   __syncthreads();
   if(threadIdx.x == 0)
   {
      CDF_BLOCK[blockIdx.x] = shared_CDF[0];
      atomicMax(&CDFblock, blockIdx.x);
   }
}

__global__ void gpu_histEq_CDFmin_02(double *CDF, double *CDF_BLOCK, 
                                     int nInpBins)
{
/* reduction method used to compute variable CDFmin
   started in kernel gpu_histEq_CDFmin_01 and
   completed in kernel gpu_histEq_CDFmin_02
*/
   // thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= 1 ) return;
   double m = 9999999999;
   for(int j = 0; j <= CDFblock; j++)
   {
      if(m > CDF_BLOCK[j]) m = CDF_BLOCK[j];
   }
   CDF_BLOCK[0] = m;
}

__global__ void gpu_histEq_LUT(double *LUT, double *CDF, double *CDF_BLOCK, 
                               int nInpBins, int nOutBins, int nPxls)
{
/* reduction method used to compute variable LUT */
   // shared variables using shared memory
   __shared__ double shared_CDF_BLOCK[1];
   // thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= nInpBins ) return;
   // thread ID for given bock size and number of blocks
   if(threadIdx.x == 0) shared_CDF_BLOCK[0] = CDF_BLOCK[0];
   __syncthreads();
   LUT[i] = ((CDF[i] - shared_CDF_BLOCK[0]) * (double)(nOutBins - 1)) / 
            ((double)nPxls - shared_CDF_BLOCK[0]);
}

__global__ void gpu_histEq_out(double *LUT, int *out, int *streamA, int nPxls)
{
/* compute array out */
   //thread ID for given bock size and number of blocks
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= nPxls ) return;
   out[i] = LUT[(int)streamA[i]];
}

extern "C" int
gpu_hist_eq(int *data, int *out, int nRows, int nCols, int nBpp, int nInpBpp, 
            int nOutBpp)
{
   int nPxls = nRows * nCols;
   int nBins = (1 << 16);
   int nOutBins = (1 << nOutBpp);
   int nInpBins = (1 << nInpBpp);
   // define and allocate device arrays
   int size_data = nRows * nCols * sizeof(int);
   int size_h = nBins * sizeof(int);
   int size_CDF = nInpBins * sizeof(int);
   int size_LUT = nInpBins * sizeof(int);
   int size_out = nPxls * sizeof(int);
   int *dev_streamA;
   int *dev_h;
   double *dev_CDF;
   double *dev_LUT;
   int *dev_out;
   double *dev_CDF_BLOCK;
   int size_CDF_BLOCK = (nInpBins / BLOCK_SIZE + 1) * sizeof(double);
   CUDA_SAFE_MALLOC(dev_streamA, size_data );
   CUDA_SAFE_MALLOC(dev_h, size_h );
   CUDA_SAFE_MALLOC(dev_CDF, size_CDF );
   CUDA_SAFE_MALLOC(dev_LUT, size_LUT );
   CUDA_SAFE_MALLOC(dev_out, size_out );
   CUDA_SAFE_MALLOC(dev_CDF_BLOCK, size_CDF_BLOCK );
   // copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemcpy(dev_streamA, data, size_data, cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemset(dev_h, 0, size_h ));
   CUDA_SAFE(cudaMemset(dev_CDF, 0, size_CDF ));
   CUDA_SAFE(cudaMemset(dev_LUT, 0, size_LUT ));
   CUDA_SAFE(cudaMemset(dev_out, 0, size_out ));
   CUDA_SAFE(cudaMemset(dev_CDF_BLOCK, 0, size_CDF_BLOCK ));
   // compute number of blocks for kernel launch
   int num_blocks = ( nPxls ) / BLOCK_SIZE + 1;
/* kernel gpu_hist perform the tasks equivalent to function hist in C code */
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_hist<<<num_blocks,BLOCK_SIZE>>>(dev_streamA, dev_h, nRows, nCols, nBins);
   // Synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // compute number of blocks for kernel launch
   num_blocks = ( nInpBins ) / BLOCK_SIZE + 1;

/* combined kernels 
 gpu_histEq_CDF_01
 gpu_histEq_CDF_02
 gpu_histEq_CDFmin_01
 gpu_histEq_CDFmin_02
 gpu_histEq_LUT
 gpu_histEq_out
 perform the tasks equivalent to function histEq in C code */

   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_histEq_CDF_01<<<num_blocks,BLOCK_SIZE>>>(dev_h, dev_CDF, dev_CDF_BLOCK, 
                                                nInpBins);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // compute number of blocks for kernel launch
   num_blocks = ( nInpBins ) / BLOCK_SIZE + 1;
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_histEq_CDF_02<<<num_blocks,BLOCK_SIZE>>>(dev_CDF, dev_CDF_BLOCK, 
                                                nInpBins);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // compute number of blocks for kernel launch
   num_blocks = ( nInpBins + BLOCK_SIZE - 1) / BLOCK_SIZE;
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_histEq_CDFmin_01<<<num_blocks,BLOCK_SIZE>>>(dev_CDF, dev_CDF_BLOCK, 
                                                   nInpBins);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_histEq_CDFmin_02<<<1,1>>>(dev_CDF, dev_CDF_BLOCK, nInpBins);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // compute number of blocks for kernel launch
   num_blocks = ( nInpBins + BLOCK_SIZE - 1) / BLOCK_SIZE;
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_histEq_LUT<<<num_blocks,BLOCK_SIZE>>>(dev_LUT, dev_CDF, dev_CDF_BLOCK, 
                                             nInpBins, nOutBins, nPxls);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // compute number of blocks for kernel launch
   num_blocks = ( nPxls + BLOCK_SIZE - 1) / BLOCK_SIZE;
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_histEq_out<<<num_blocks,BLOCK_SIZE>>>(dev_LUT, dev_out, dev_streamA, 
                                             nPxls);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(out, dev_out,size_out, cudaMemcpyDeviceToHost));
   // free device variables
   cudaFree(dev_streamA);
   cudaFree(dev_h);
   cudaFree(dev_CDF);
   cudaFree(dev_LUT);
   cudaFree(dev_out);
   cudaFree(dev_CDF_BLOCK);
   return 0;
}
