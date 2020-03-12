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
	 This CUDA code is implementation of the kernel - 2d convolution
	 as per C code in file 2dconv.c provided with header as above

	 The cuda implementation done by the team at PNNL
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../inc/2d_convolution.h"

/*
	 define block size for kernel launch, could be variable
	 tested from a value of 64 upto a maximum of 1024
 */

#define BLOCK_SIZE 512
#define MAX_NFILTERROWS 10
#define MAX_NFILTERCOLS 10

// Helpful utility macros

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

// Array filter stored in constant memory
__constant__ float filter[MAX_NFILTERROWS*MAX_NFILTERCOLS]; // 800 bytes; 0.78125 KB

__global__ void compute_conv(int nRows, int nCols, int nFilterRows, 
                             int nFilterCols, int rowOffset, int colOffset, 
                             int rowBegIndex, int colBegIndex, int *tmpBuf, 
                             int *out, float normFactor)
{
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x * blockDim.x) + threadIdx.x;
   if ( thread_number >= (nRows + rowOffset - rowBegIndex) * 
                         (nCols + colOffset - colBegIndex) ) 
      return;
   int row_ = thread_number / (nCols + colOffset - colBegIndex);
   int col_ = thread_number - row_ * (nCols + colOffset - colBegIndex); 
   int row = row_ + rowBegIndex;
   int col = col_ + colBegIndex;
   if ( row >= nRows + rowOffset) return;
   if (col >= nCols + colOffset ) return;
   float sum = 0.0;
   int m, n;
   int i, j;
   int pxlPos, fltPos;
   sum = 0;
   m = 0;
#pragma unroll 
   for (i = row - rowOffset; i <= row + rowOffset; i++)
   {
      n = 0;
      for (j = col - colOffset; j <= col + colOffset; j++)
      {
         pxlPos = i * (nCols + nFilterCols) + j;
         fltPos = m * nFilterCols + n;
         sum += ( (float)tmpBuf[pxlPos] * filter[fltPos]);
         n++;
      }
      m++;
   }
   out[(row - rowBegIndex) * nCols + (col - colBegIndex)] = 
      (int)(sum / normFactor);
}

extern "C" int
gpu_conv2d (int *in, int *out, int nRows, int nCols, float *_filter, 
            double normFactor, int nFilterRows, int nFilterCols)
{
   int row = 0;
   int rowOffset = nFilterRows >> 1;
   int colOffset = nFilterCols >> 1;
   int rowBegIndex = rowOffset;
   int colBegIndex = colOffset;
   // Get the size of float array filter and copy to location in constant memory
   int size_filter = nFilterRows * nFilterCols * sizeof(int);
   CUDA_SAFE(cudaMemcpyToSymbol(filter, _filter, size_filter, 0, 
             cudaMemcpyHostToDevice));
   int size_in_out = ((nRows + nFilterRows) * (nCols + nFilterCols)*sizeof(int) );
   int size_out = (nRows * nCols *sizeof(int) );
   // Allocate memory for variable tmpBuf on device and initialize to zero
   int *dev_tmpBuf;
   CUDA_SAFE_MALLOC(dev_tmpBuf, size_in_out );
   CUDA_SAFE(cudaMemset(dev_tmpBuf, 0, size_in_out));
   // Copy variable tmpBuf from Host to Device
   for (row = 0; row < nRows; row++)
   {
      CUDA_SAFE(cudaMemcpy(dev_tmpBuf+ (row + rowOffset) * (nCols + nFilterCols) + 
                colOffset, in+ row * nCols, nCols * sizeof(int), 
                cudaMemcpyHostToDevice));  
   }
   // Allocate memory for variable out on device and initialize to zero
   int *dev_out;
   CUDA_SAFE_MALLOC(dev_out, size_out);
   CUDA_SAFE(cudaMemset(dev_out, 0, size_out));
   // Compute number of blocks for kernel launch
   int num_blocks = (nRows + rowOffset - rowBegIndex) * (nCols + colOffset - 
                     colBegIndex) / BLOCK_SIZE + 1;
   // Launch kernel for num_blocks and BLOCK_SIZE
   compute_conv<<<num_blocks,BLOCK_SIZE>>>( nRows, nCols, nFilterRows, 
      nFilterCols, rowOffset, colOffset, rowBegIndex, colBegIndex, 
      dev_tmpBuf, dev_out, normFactor);
   // Synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // Copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(out, dev_out,size_out, cudaMemcpyDeviceToHost));
   // Free device variables
   cudaFree(dev_tmpBuf);
   cudaFree(dev_out);
   return 0;
}
