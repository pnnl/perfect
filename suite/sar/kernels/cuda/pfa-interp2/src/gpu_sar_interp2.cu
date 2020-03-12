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

#include "sar_interp2.h"
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
#define BLOCK_SIZE 128

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) {  \
   printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}

#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))
#define PFA_N_TSINC_POINTS_PER_SIDE ((T_PFA - 1)/2)

__global__ void gpu_sar_interp2_kernel_011 (double *input_coords, float *input_spacing_avg)
{
   // thread ID for given bock size and number of blocks
   int r = (blockIdx.x * blockDim.x) + threadIdx.x; 
   if(r >= PFA_NOUT_RANGE) return;
   int p; 
   float _input_spacing_avg;
   _input_spacing_avg = 0.0f;
   for (p = 0; p < N_PULSES-1; ++p)
   {
      _input_spacing_avg += fabs(input_coords[r*N_PULSES+p+1] - 
                                 input_coords[r*N_PULSES+p]);
   }
   input_spacing_avg[r] = _input_spacing_avg;
}

__global__ void gpu_sar_interp2_kernel_02 (double *input_coords, 
   float *input_spacing_avg, double *output_coords, 
   complex *resampled, float *window, complex *data)
{
   __shared__ float shared_input_spacing_avg;
   __shared__ float shared_input_spacing_avg_inv;
   __shared__ float shared_scale_factor;

   // thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(thread_number >= PFA_NOUT_RANGE * PFA_NOUT_AZIMUTH ) return;

   int r = thread_number / PFA_NOUT_AZIMUTH;
   if(r >= PFA_NOUT_RANGE ) return;

   int p = thread_number - r * PFA_NOUT_AZIMUTH;
   if(p >= PFA_NOUT_AZIMUTH) return;

   /* shared variables are accessed only by thread-zero in order 
      to reduce memory transctions */
   if(threadIdx.x == 0) 
   {
      shared_input_spacing_avg = input_spacing_avg[r];
      shared_input_spacing_avg /= (N_PULSES-1);
      shared_input_spacing_avg_inv = 1.0f / shared_input_spacing_avg;
      shared_scale_factor = fabs(output_coords[1] - output_coords[0]) * 
                            shared_input_spacing_avg_inv;
   }
   __syncthreads();
   const double out_coord = output_coords[p];
   /* Begin part of code to find nearest */
   int nearest;
   int left_ind, right_ind, mid_ind;
   double left_val, right_val, mid_val;

   /*
    * We assume for simplicity that the input coordinates are
    * monotonically increasing.
    */
   // assert(PFA_NOUT_RANGE > 1 && input_coords[1] > input_coords[0]);
   if(PFA_NOUT_RANGE <= 1 || input_coords[1] <= input_coords[0])
   {
      printf("ASSERT FAIL for PFA_NOUT_RANGE <= 1 "); 
      printf("or input_coords[1] <= input_coords[0] \n");
      return;
   }
   left_ind = 0;
   right_ind = PFA_NOUT_RANGE-1;
   mid_ind = (left_ind+right_ind)/2;
   left_val = input_coords[r*N_PULSES+left_ind];
   right_val = input_coords[r*N_PULSES+right_ind];
   mid_val = input_coords[r*N_PULSES+mid_ind];

   // target_coord becomes output_coords[p]
   // input_coord need to have 2-dimensions

   if (output_coords[p] < left_val || output_coords[p] > right_val)
   {
      nearest = -1;
   }
   else
   {
      while (right_ind - left_ind > 1)
      {
         if (output_coords[p] <= mid_val)
         {
            right_ind = mid_ind;
            right_val = mid_val;
         }
         else
         {
            left_ind = mid_ind;
            left_val = mid_val;
         }
         mid_ind = (left_ind+right_ind)/2;
         mid_val = input_coords[r*N_PULSES+mid_ind]; 
      }
      nearest = mid_ind;
   }

   /* End part of code to find nearest */

   if (nearest < 0)
   {
      resampled[p*PFA_NOUT_RANGE+r].re = 0.0f;
      resampled[p*PFA_NOUT_RANGE+r].im = 0.0f;
   }
   else
   {
      /* find_nearest_azimuth_coord should never return a value >= N_PULSES */
      if(nearest >= N_PULSES )
      {
         printf("ASSERT FAIL for nearest < N_PULSES \n");
         return;
      }
      /*
       * out_coord is bounded in [nearest, nearest+1], so we check
       * which of the two input coordinates is closest.
       */
      if (fabs(out_coord-input_coords[r*N_PULSES+nearest+1]) < 
          fabs(out_coord-input_coords[r*N_PULSES+nearest]))
      {
         nearest = nearest + 1;
      }

      int k, pmin, pmax, window_offset;
      complex accum;
      float sinc_arg, sinc_val;

      pmin = nearest - PFA_N_TSINC_POINTS_PER_SIDE;
      if (pmin < 0) { pmin = 0; }
      pmax = nearest + PFA_N_TSINC_POINTS_PER_SIDE;
      if (pmax >= N_PULSES) { pmax = N_PULSES-1; }
      window_offset = 0;
      if (nearest - PFA_N_TSINC_POINTS_PER_SIDE < 0)
      {
         window_offset = PFA_N_TSINC_POINTS_PER_SIDE - nearest;
      }
      accum.re = accum.im = 0.0f;
      for (k = pmin; k <= pmax; ++k)
      {
         sinc_arg = (out_coord - input_coords[r*N_PULSES+k]) * 
                     shared_input_spacing_avg_inv;
         if(sinc_arg == 0)
         {
            sinc_val = 1.0f;
         }
         else
         {
            const float arg = M_PI * sinc_arg;
            sinc_val = (float) sin(arg) / arg;
         }
         accum.re += sinc_val * window[window_offset+(k-pmin)] * 
                     data[k*PFA_NOUT_RANGE+r].re;
         accum.im += sinc_val * window[window_offset+(k-pmin)] * 
                     data[k*PFA_NOUT_RANGE+r].im;
      }
      resampled[p*PFA_NOUT_RANGE+r].re = shared_scale_factor * accum.re;
      resampled[p*PFA_NOUT_RANGE+r].im = shared_scale_factor * accum.im;

   }
}

extern "C" int gpu_sar_interp2(
   complex (* resampled)[PFA_NOUT_RANGE],
   complex (* const data)[PFA_NOUT_RANGE],
   const float *window,
   double (* const input_coords)[N_PULSES],
   const double *output_coords)
{
   double *input_coords_01 = &input_coords[0][0];
   // define and allocate device arrays
   float *dev_input_spacing_avg;
   double *dev_input_coords;
   double *dev_output_coords;
   complex *dev_resampled;
   float *dev_window;
   complex *dev_data;

   const size_t num_data_elements = N_PULSES * PFA_NOUT_RANGE;
   const size_t num_resampled_elements = PFA_NOUT_AZIMUTH * PFA_NOUT_RANGE;
   const size_t num_window_elements = T_PFA;

   CUDA_SAFE_MALLOC(dev_input_spacing_avg, PFA_NOUT_RANGE * sizeof(float) );
   CUDA_SAFE_MALLOC(dev_input_coords, num_data_elements * sizeof(double) );
   CUDA_SAFE_MALLOC(dev_output_coords, PFA_NOUT_RANGE * sizeof(double) );

   CUDA_SAFE(cudaMalloc((void **)&dev_resampled, 
      num_resampled_elements * sizeof(complex)));
   CUDA_SAFE_MALLOC(dev_window, num_window_elements * sizeof(float));
   CUDA_SAFE(cudaMalloc((void **)&dev_data, 
      num_data_elements * sizeof(complex)));
   
   // copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemset(dev_input_spacing_avg, 0, 
      PFA_NOUT_RANGE * sizeof(float) ));
   CUDA_SAFE(cudaMemcpy(dev_input_coords, input_coords_01,  
      num_data_elements * sizeof(double), cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_output_coords, output_coords,  
      PFA_NOUT_RANGE * sizeof(double), cudaMemcpyHostToDevice));

   CUDA_SAFE(cudaMemset(dev_resampled, 0, 
      num_resampled_elements * sizeof(complex)));
   CUDA_SAFE(cudaMemcpy(dev_window, window, 
      num_window_elements * sizeof(float), cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_data, data, num_data_elements * sizeof(complex), 
      cudaMemcpyHostToDevice));

   if (N_PULSES == 0 || N_RANGE == 0 || PFA_NOUT_RANGE == 0)
   {
      return 0;
   }
   assert(N_RANGE > 1 && PFA_NOUT_RANGE > 1);
   if(N_RANGE <= 1 || PFA_NOUT_RANGE <= 1){
      printf("ASSERT FAILED FOR N_RANGE > 1 or PRA_NOUT_RANGE > 1 \n");
      return 0;
   }
   int num_blocks;
   // compute number of blocks for kernel launch
   num_blocks = ( PFA_NOUT_RANGE + BLOCK_SIZE - 1) / BLOCK_SIZE;
   /* kernel gpu_sar_interp2_kernel_011 computes the values of 
      input_spacing_avg */

   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_sar_interp2_kernel_011<<<num_blocks,BLOCK_SIZE>>>(dev_input_coords, 
      dev_input_spacing_avg );

   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());

   // compute number of blocks for kernel launch
   num_blocks = ( PFA_NOUT_RANGE * PFA_NOUT_AZIMUTH + BLOCK_SIZE - 1) / 
                BLOCK_SIZE;

   /* kernel gpu_sar_interp2_kernel_02 computes the values of output coord */

   // launch kernel for num_blocks and BLOCK_SIZE
   gpu_sar_interp2_kernel_02<<<num_blocks,BLOCK_SIZE>>>(dev_input_coords, 
      dev_input_spacing_avg, dev_output_coords, 
      dev_resampled, dev_window, dev_data);

   // synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());

   // copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(resampled, dev_resampled, 
            (num_resampled_elements * sizeof(complex)), 
             cudaMemcpyDeviceToHost));

   // free device variables
   CUDA_SAFE(cudaFree(dev_input_coords));
   CUDA_SAFE(cudaFree(dev_output_coords));
   CUDA_SAFE(cudaFree(dev_resampled));
   CUDA_SAFE(cudaFree(dev_window));
   CUDA_SAFE(cudaFree(dev_data));
   return 0;
}
