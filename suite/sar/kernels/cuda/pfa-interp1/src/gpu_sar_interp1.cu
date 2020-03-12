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


#include "sar_interp1.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "../lib/sar_utils.h"

/*
   define block size for kernel launch, could be variable
   tested from a value of 128 upto a maximum of 1024
   based on a different grid size, block size could be changed
 */
#define BLOCK_SIZE 256

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

__global__ void gpu_sar_interp1_kernel (double *input_coords_start, 
   double *input_coords_spacing, double *output_coords, 
   complex *resampled, float *window, complex *data)
{
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(thread_number >= N_PULSES * PFA_NOUT_RANGE ) return;
   int p = (thread_number) / PFA_NOUT_RANGE;
   if(p >= N_PULSES) return;
   int r = thread_number - p * PFA_NOUT_RANGE;
   if(r >= PFA_NOUT_RANGE ) return;
   double input_spacing, input_start, input_spacing_inv;
   float scale_factor;
   input_start = input_coords_start[p];
   input_spacing = input_coords_spacing[p];
   input_spacing_inv = 1.0 / input_spacing;
   scale_factor = fabs(output_coords[1] - output_coords[0]) * input_spacing_inv;
   const double out_coord = output_coords[r];
   int nearest;
   double x;
   int k, rmin, rmax, window_offset;
   const int PFA_N_TSINC_POINTS_PER_SIDE = (T_PFA - 1)/2;
   complex accum;
   float sinc_arg, sinc_val, win_val;
   //find nearest range coordinate
   /*
    * Test for the target coordinate being out-of-bounds with respect to
    * the input coordinates.
    */
   if (output_coords[r] < input_start || output_coords[r] >= 
      (input_start + (N_RANGE-1)*input_spacing))
   {
      nearest = -1;
   }
   else
   {
      x = (output_coords[r] - input_start) * input_spacing_inv;
      nearest = (int)(x + 0.5);
   }
   if (nearest < 0)
   {
      resampled[p*PFA_NOUT_RANGE+r].re = 0.0f;
      resampled[p*PFA_NOUT_RANGE+r].im = 0.0f;
   }
   else
   {
      /* find_nearest_range_coord should never return a value >= N_RANGE */
      if(nearest >= N_RANGE){
         printf("ASSERT FAILED FOR nearest < N_RANGE\n");
         return;
      }
      /*
       * out_coord is bounded in [nearest, nearest+1], so we check
       * which of the two input coordinates is closest.
       */
      if (fabs(out_coord - (input_start + (nearest+1)*input_spacing)) <
          fabs(out_coord - (input_start + (nearest)*input_spacing)))
      {
         nearest = nearest + 1;
      }	
      rmin = nearest - PFA_N_TSINC_POINTS_PER_SIDE;
      if (rmin < 0) { rmin = 0; }
      rmax = nearest + PFA_N_TSINC_POINTS_PER_SIDE;
      if (rmax >= N_RANGE) { rmax = N_RANGE-1; }
      window_offset = 0;
      if (nearest - PFA_N_TSINC_POINTS_PER_SIDE < 0)
      {
         window_offset = PFA_N_TSINC_POINTS_PER_SIDE - nearest;
      }
      accum.re = accum.im = 0.0f;
      for (k = rmin; k <= rmax; ++k)
      {
         win_val = window[window_offset+(k-rmin)];
         sinc_arg = (out_coord - (input_start+k*input_spacing)) * 
                     input_spacing_inv;
         if(sinc_arg == 0)  
         {
            sinc_val = 1.0f;
         }
         else
         {
            const float arg = M_PI * sinc_arg;
            sinc_val = (float) sin(arg) / arg;
         }
         accum.re += sinc_val * win_val * data[p*N_RANGE+k].re;
         accum.im += sinc_val * win_val * data[p*N_RANGE+k].im;
      }
      resampled[p*PFA_NOUT_RANGE+r].re = scale_factor * accum.re;
      resampled[p*PFA_NOUT_RANGE+r].im = scale_factor * accum.im;
   }
}

extern "C" int gpu_sar_interp1(
   complex (* resampled)[PFA_NOUT_RANGE],
   complex (* const data)[N_RANGE],
   const float *window,
   const double input_coords_start[N_PULSES],
   const double input_coords_spacing[N_PULSES],
   const double output_coords[PFA_NOUT_RANGE])
{
   if (N_PULSES == 0 || N_RANGE == 0 || PFA_NOUT_RANGE == 0)
   {
      return 0;
   }		
   assert(N_RANGE > 1 && PFA_NOUT_RANGE > 1);
   if(N_RANGE <= 1 || PFA_NOUT_RANGE <= 1){
      printf("ASSERT FAILED FOR N_RANGE > 1 or PRA_NOUT_RANGE > 1 \n");
      return 0;
   }
   complex *data_01 = &data[0][0];
   // define and allocate device arrays
   double *dev_input_coords_start;
   double *dev_input_coords_spacing;
   double *dev_output_coords;
		
   complex *dev_resampled;
   float *dev_window;
   complex *dev_data;

   const size_t num_data_elements = N_PULSES * N_RANGE;
   const size_t num_resampled_elements = N_PULSES * PFA_NOUT_RANGE;
   const size_t num_window_elements = T_PFA;

   CUDA_SAFE_MALLOC(dev_input_coords_start, N_PULSES * sizeof(double));
   CUDA_SAFE_MALLOC(dev_input_coords_spacing, N_PULSES * sizeof(double));
   CUDA_SAFE_MALLOC(dev_output_coords, PFA_NOUT_RANGE * sizeof(double));

   CUDA_SAFE(cudaMalloc((void **)&dev_resampled, 
      num_resampled_elements * sizeof(complex)));
   CUDA_SAFE_MALLOC(dev_window, num_window_elements * sizeof(float));
   CUDA_SAFE(cudaMalloc((void **)&dev_data, 
      num_data_elements * sizeof(complex)));

   //copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemcpy(dev_input_coords_start, input_coords_start, 
      N_PULSES * sizeof(double), cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_input_coords_spacing, input_coords_spacing, 
      N_PULSES * sizeof(double), cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_output_coords, output_coords, 
      PFA_NOUT_RANGE * sizeof(double), cudaMemcpyHostToDevice));

   CUDA_SAFE(cudaMemset(dev_resampled, 0, 
      num_resampled_elements * sizeof(complex)));

   CUDA_SAFE(cudaMemcpy(dev_window, window, 
      num_window_elements * sizeof(float), cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_data, data_01, 
      num_data_elements * sizeof(complex), cudaMemcpyHostToDevice));

   int num_blocks;

   // compute number of blocks for kernel launch
   num_blocks = ( N_PULSES * PFA_NOUT_RANGE + BLOCK_SIZE - 1) / BLOCK_SIZE;

   /* kernel gpu_sar_interp1_kernel performs the tasks equivalent to 
      function sar_interp1 */

   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_sar_interp1_kernel<<<num_blocks,BLOCK_SIZE>>>(dev_input_coords_start, 
      dev_input_coords_spacing, dev_output_coords, 
      dev_resampled, dev_window, dev_data);

   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());

   // copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(resampled, dev_resampled, 
      (num_resampled_elements * sizeof(complex)), cudaMemcpyDeviceToHost));

   // free device variables
   CUDA_SAFE(cudaFree(dev_input_coords_start));
   CUDA_SAFE(cudaFree(dev_input_coords_spacing));
   CUDA_SAFE(cudaFree(dev_output_coords));

   CUDA_SAFE(cudaFree(dev_resampled));
   CUDA_SAFE(cudaFree(dev_window));
   CUDA_SAFE(cudaFree(dev_data));

   return 0;
}
