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
   This CUDA code is implementation of the kernel - discrete wavelet transform
   as per C code in file dwt.c provided with header as above

   The cuda implementation done by the team at PNNL
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../inc/dwt53.h"

/*
   define block size for kernel launch, could be variable
   tested from a value of 64 upto a maximum of 1024
 */

#define BLOCK_SIZE 64 
#define MAX_NFILTERROWS 10
#define MAX_NFILTERCOLS 10

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) \
        { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

__global__ void dwt53_transpose_01(int nrows, int ncols, int *data, int *data2)
{
   // thread ID for given bock size and number of blocks
   // integer variables i and j have meaning corresponding to the C code
   int i = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( i >= nrows ) return;
   int j;
   int cur;
   /* Predict the odd pixels using linear interpolation of the even pixels */
   for (j = 1; j < ncols - 1; j += 2)
   {
      cur = i * ncols + j;
#ifdef USE_SHIFT
      data[cur] -= (data[cur - 1] + data[cur + 1]) >> 1;
#else
      data[cur] -= (int)(0.5 * (data[cur - 1] + data[cur + 1]));
#endif
   }
/* The last odd pixel only has its left neighboring even pixel */
   cur = i * ncols + ncols - 1;
   data[cur] -= data[cur - 1];
/* Update the even pixels using the odd pixels
 * to preserve the mean value of the pixels
 */
   for (j = 2; j < ncols; j += 2)
   {
      cur = i * ncols + j;
#ifdef USE_SHIFT
      data[cur] += (data[cur - 1] + data[cur + 1]) >> 2;
#else
      data[cur] += (int)(0.25 * (data[cur - 1] + data[cur + 1]));
#endif
   }
/* The first even pixel only has its right neighboring odd pixel */
   cur = i * ncols + 0;
#ifdef USE_SHIFT
   data[cur] += data[cur + 1] >> 1;
#else
   data[cur] += (int)(0.5 * data[cur + 1]);
#endif
   for (j = 0; j < ncols / 2; j++)
   {
      data2[j * nrows + i] = data[i * ncols + 2 * j];
      data2[(j + ncols / 2)* nrows + i] = data[i * ncols + 2 * j + 1];
   }
}

extern "C" int
gpu_dwt53(int *data, int nrows, int ncols)
{
   // define and allocate device arrays
   int size_data = nrows * ncols * sizeof(int);
   int *dev_data;
   int *dev_data2;
   CUDA_SAFE_MALLOC(dev_data, size_data );
   CUDA_SAFE_MALLOC(dev_data2, size_data );

   // copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemcpy(dev_data, data, size_data, cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemset(dev_data2, 0, size_data));

/* First do all rows; This function will transpose the data 
 * as it performs its final shuffling
 */
   // compute number of blocks for kernel launch
   int num_blocks = ( nrows + BLOCK_SIZE -1 ) / BLOCK_SIZE;
   // launch kernel for num_blocks and BLOCK_SIZE
   dwt53_transpose_01<<<num_blocks,BLOCK_SIZE>>>(nrows, ncols, dev_data, 
                                                 dev_data2);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
/* We next do all the columns (they are now the rows) */
   // compute number of blocks for kernel launch
   num_blocks = ( ncols ) / BLOCK_SIZE + 1;
   // launch kernel for num_blocks and BLOCK_SIZE
   dwt53_transpose_01<<<num_blocks,BLOCK_SIZE>>>(ncols, nrows, dev_data2, 
                                                 dev_data);
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(data, dev_data,size_data, cudaMemcpyDeviceToHost));
   // free device variables
   cudaFree(dev_data);
   cudaFree(dev_data2);
   return 0;
}
