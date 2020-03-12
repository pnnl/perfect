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

#include "sar_backprojection.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/*
   define block size for kernel launch, could be variable
   tested from a value of 128 upto a maximum of 512
   based on a different grid size, block size could be changed
 */

#define BLOCK_SIZE 128

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) {                    \
                     printf("CUDA CALL FAILED AT %d\n", __LINE__ ); \
                     exit(1);}                                      
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

__global__ void gpu_backprojection_kernel (complex *image, complex *data, 
           double * platpos, double ku, double R0, double dR, double dxdy, 
           double z0)
   {
   // thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(thread_number >= BP_NPIX_Y * BP_NPIX_X ) return;
   // integer variable iy is equivalent to variable iy in C code
   int iy = (thread_number) / BP_NPIX_X;
   if(iy >= BP_NPIX_Y) return;
   // integer variable ix is equivalent to variable ix in C code
   int ix = thread_number - iy * BP_NPIX_X;
   if(ix >= BP_NPIX_X ) return;
   int p;
   const double dR_inv = 1.0/dR;
   complex accum;
   const double py = (-BP_NPIX_Y/2.0 + 0.5 + iy) * dxdy;
   const double px = (-BP_NPIX_X/2.0 + 0.5 + ix) * dxdy;
   accum.re = accum.im = 0.0f;
   for (p = 0; p < N_PULSES; ++p)
   {
      /* calculate the range R from the platform to this pixel */
      const double xdiff = platpos[p*3+0] - px;
      const double ydiff = platpos[p*3+1] - py;
      const double zdiff = platpos[p*3+2] - z0;
      const double R = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);
      /* convert to a range bin index */
      const double bin = (R-R0)*dR_inv;
      if (bin >= 0 && bin <= N_RANGE_UPSAMPLED-2)
      {
         complex sample, matched_filter, prod;
         /* interpolation range is [bin_floor, bin_floor+1] */
         const int bin_floor = (int) bin;
         /* interpolation weight */
         const float w = (float) (bin - (double) bin_floor);
         /* linearly interpolate to obtain a sample at bin */
         sample.re = (1.0f-w)*data[p*N_RANGE_UPSAMPLED+bin_floor].re + 
                     w*data[p*N_RANGE_UPSAMPLED+bin_floor+1].re;
         sample.im = (1.0f-w)*data[p*N_RANGE_UPSAMPLED+bin_floor].im + 
                     w*data[p*N_RANGE_UPSAMPLED+bin_floor+1].im;
         /* compute the complex exponential for the matched filter */
         matched_filter.re = cos(2.0 * ku * R);
         matched_filter.im = sin(2.0 * ku * R);
         /* scale the interpolated sample by the matched filter */
         prod.re = sample.re * matched_filter.re - 
                   sample.im * matched_filter.im;
         prod.im = sample.re * matched_filter.im + 
                   sample.im * matched_filter.re;
         /* accumulate this pulse's contribution into the pixel */
         accum.re += prod.re;
         accum.im += prod.im;
      }
   }
   image[iy*BP_NPIX_X+ix] = accum;
}

extern "C" int gpu_sar_backprojection(
   complex image[BP_NPIX_Y][BP_NPIX_X],
   complex (* const data)[N_RANGE_UPSAMPLED],
   const position platpos[N_PULSES],
   double ku,
   double R0,
   double dR,
   double dxdy,
   double z0)
{
   complex *data_01 = &data[0][0];
   // define and allocate device arrays
   complex *dev_image;
   complex *dev_data;
   double *dev_platpos;
   const size_t num_image_elements = BP_NPIX_Y * BP_NPIX_X;
   const size_t num_data_elements = N_PULSES * N_RANGE_UPSAMPLED;
   CUDA_SAFE(cudaMalloc((void **)&dev_image, num_image_elements * 
             sizeof(complex)));
   CUDA_SAFE(cudaMalloc((void **)&dev_data, num_data_elements * 
             sizeof(complex)));
   CUDA_SAFE_MALLOC(dev_platpos, N_PULSES * sizeof(position));
   // copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemset(dev_image, 0, num_image_elements * sizeof(complex) ));
   CUDA_SAFE(cudaMemcpy(dev_data, data_01, num_data_elements * sizeof(complex),
                        cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_platpos, platpos, N_PULSES * sizeof(position), 
                        cudaMemcpyHostToDevice));
   // compute number of blocks for kernel launch
   int num_blocks = ( BP_NPIX_Y * BP_NPIX_X + BLOCK_SIZE - 1) / BLOCK_SIZE;
   /* kernel gpu_backprojection_kernel perform the tasks equivalent to function
      sar_backprojection in C code */
   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_backprojection_kernel<<<num_blocks,BLOCK_SIZE>>>(dev_image, dev_data, 
                            dev_platpos, ku, R0, dR, dxdy, z0 );
   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   // copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(image, dev_image, (num_image_elements * 
             sizeof(complex)), cudaMemcpyDeviceToHost));
   // free device variables
   CUDA_SAFE(cudaFree(dev_image));
   CUDA_SAFE(cudaFree(dev_data));
   CUDA_SAFE(cudaFree(dev_platpos));
   return 0;
}
