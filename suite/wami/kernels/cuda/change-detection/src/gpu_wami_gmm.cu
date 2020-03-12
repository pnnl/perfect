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

#include "wami_utils.h"
#include <assert.h>

#define BLOCK_SIZE 32

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

#define ONE_OVER_SQRT_TWO_PI (1.0f / sqrt(2.0f * M_PI))

__global__ void gpu_wami_gmm_kernel(unsigned char* foreground, float *mu, 
   float *sigma, float *weight, unsigned short *frame)
{
   __shared__ float shared_mu[BLOCK_SIZE * 5];
   __shared__ float shared_sigma[BLOCK_SIZE * 5];
   __shared__ float shared_weight[BLOCK_SIZE * 5];

   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * 
      blockDim.y + (threadIdx.y * blockDim.x) + threadIdx.x;

   if(thread_number >= WAMI_GMM_IMG_NUM_ROWS * WAMI_GMM_IMG_NUM_COLS) return;

   //thread ID for given bock size and number of blocks
   int col = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(col >= WAMI_GMM_IMG_NUM_COLS) return;

   int row = blockIdx.y;
   if(row >= WAMI_GMM_IMG_NUM_ROWS) return;

   int col_0 = (blockIdx.x * blockDim.x);
   int k;

   //copy block of most used variable arrays to shared memory
   for(k = 0; k < WAMI_GMM_NUM_MODELS; k++)
   {
      shared_mu[k*blockDim.x + threadIdx.x] = mu[row*(WAMI_GMM_IMG_NUM_COLS*
         WAMI_GMM_NUM_MODELS) + col_0*WAMI_GMM_NUM_MODELS + k*blockDim.x + 
         threadIdx.x];
      shared_sigma[k*blockDim.x + threadIdx.x] = sigma[row*
         (WAMI_GMM_IMG_NUM_COLS*WAMI_GMM_NUM_MODELS) + col_0*WAMI_GMM_NUM_MODELS
         + k*blockDim.x + threadIdx.x];
      shared_weight[k*blockDim.x + threadIdx.x] = weight[row*
         (WAMI_GMM_IMG_NUM_COLS*WAMI_GMM_NUM_MODELS) + col_0*WAMI_GMM_NUM_MODELS
         + k*blockDim.x + threadIdx.x];
   }
   __syncthreads();
   const float STDEV_THRESH = 2.5f;
   const float INIT_STDEV = 80.0f;
   const float alpha = 0.01f; /* Learning rate */
   const float INIT_WEIGHT = 0.01f;
   const float BACKGROUND_THRESH = 0.9f;
   const unsigned short pixel = frame[row*WAMI_GMM_IMG_NUM_COLS + col];
   int match = -1;
   float sum = 0.0f, norm = 0.0f;
   int sorted_position = 0;
   for (k = 0; k < WAMI_GMM_NUM_MODELS; ++k)
   {
   /*
    * C89 does not include fabsf(), so using the double-precision
    * fabs() function will unnecessarily type-convert to double.
    */
      if (fabs(pixel - shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + k])/
          shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + k] < STDEV_THRESH)
      {
         match = k;
         break;
      }
   }
   /* Update the weights for all models */
   for (k = 0; k < WAMI_GMM_NUM_MODELS; ++k)
   {
      if (k == match)
      {
      /* A model matched, so update its corresponding weight. */
         shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + match] += alpha * 
            (1.0f - shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + match]);
      }
      else
      {
      /* Non-matching models have their weights reduced */
         shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + k] *= (1.0f - alpha);
      }
   }

   if (match < 0)
   {
   /*
    * No distribution matched; replace the least likely distribution.
    * We keep the entries sorted by significance, so the last entry
    * is also the least likely.  We do this after updating weights
    * above so that the initial weight is not immediately down-weighted,
    * although that means that one update above was wasted. That
    * update could be avoided.
    */
      shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + (WAMI_GMM_NUM_MODELS-1)] = 
         (float) pixel;
      shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + (WAMI_GMM_NUM_MODELS-1)] = 
         INIT_STDEV;
      shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + (WAMI_GMM_NUM_MODELS-1)] =
         INIT_WEIGHT;
   }
   /* Normalize weights */
   for (k = 0; k < WAMI_GMM_NUM_MODELS; ++k)
   {
      sum += shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + k];
   }
   assert(sum != 0.0f);
   norm = 1.0f / sum;
   for (k = 0; k < WAMI_GMM_NUM_MODELS; ++k)
   {
      shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + k] *= norm;
   }
   /* Update mu and sigma for the matched distribution, if any */
   if (match >= 0)
   {
      const float mu_k = shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + match];
      const float sigma_k = shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + 
         match];
      const float sigma_k_inv = 1.0f / sigma_k;
   /*
    * C89 does not include a single-precision expf() exponential function,
    * so we instead use the double-precision variant exp(). Same for sqrt()
    * below.
    */
      const float rho = alpha * (ONE_OVER_SQRT_TWO_PI * sigma_k_inv) *
         exp( -1.0f * (pixel-mu_k)*(pixel-mu_k) / (2.0f * sigma_k * sigma_k) );
      shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + match] = (1.0f - rho) * 
         mu_k + rho * pixel;
      shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + match] = sqrt(
         (1.0f - rho) * sigma_k * sigma_k +
         rho * (pixel-shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + match]) * 
         (pixel-shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + match]));
      assert(shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + match] > 0);
   }
/*
 * weight and sigma for the matched (or new) distribution are the only
 * values that may have changed, so we find the correct location of that
 * new value in the sorted list.  Matches lead to more evidence, and thus
 * higher weight and lower sigma, so we only need to sort "higher".
 */
   sorted_position = 0;
   if (match != 0)
   {
      const int sort_from = (match >= 0) ? match : WAMI_GMM_NUM_MODELS-1;
      const float new_significance = shared_weight[threadIdx.x*
         WAMI_GMM_NUM_MODELS + sort_from] / 
         shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + sort_from];
      float other_significance, new_mu, new_sigma, new_weight;
      for (k = sort_from-1; k >= 0; --k)
      {
         other_significance = shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + k]
            / shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + k];
         if (new_significance <= other_significance)
         {  
            break;
         }
      }
      if (k == 0)
      {
         if (other_significance >= new_significance)
         {
            sorted_position = 1;
         }
         else
         {
            sorted_position = 0;
         }
      }
      else
      {
         sorted_position = k + 1;
      }
      new_mu = shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + sort_from];
      new_sigma = shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + sort_from];
      new_weight = shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + sort_from];
      for (k = sort_from; k > sorted_position; --k)
      {
         shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + k] = shared_mu[threadIdx.x*
            WAMI_GMM_NUM_MODELS + (k-1)];
         shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + k] = shared_sigma[
            threadIdx.x*WAMI_GMM_NUM_MODELS + (k-1)];
         shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + k] = shared_weight[
            threadIdx.x*WAMI_GMM_NUM_MODELS + (k-1)];
      }
      shared_mu[threadIdx.x*WAMI_GMM_NUM_MODELS + sorted_position] = new_mu;
      shared_sigma[threadIdx.x*WAMI_GMM_NUM_MODELS + sorted_position] = 
         new_sigma;
      shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + sorted_position] = 
         new_weight;
   }

   /* Now, we need to determine if this pixel is foreground or background. */
   float cumsum = shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + 0];
   int B = 0;
   while (B < WAMI_GMM_NUM_MODELS-1 && cumsum <= BACKGROUND_THRESH)
   {
      cumsum += shared_weight[threadIdx.x*WAMI_GMM_NUM_MODELS + (++B)];
   }
   foreground[row*(WAMI_GMM_IMG_NUM_COLS) + col] = (sorted_position > B);
   __syncthreads();
   //copy block of most used variable arrays back from shared memory
   for(k = 0; k < WAMI_GMM_NUM_MODELS; k++)
   {
      mu[row*(WAMI_GMM_IMG_NUM_COLS*WAMI_GMM_NUM_MODELS) + 
         col_0*WAMI_GMM_NUM_MODELS + k*blockDim.x + threadIdx.x] = shared_mu[
         k*blockDim.x + threadIdx.x];
      sigma[row*(WAMI_GMM_IMG_NUM_COLS*WAMI_GMM_NUM_MODELS) + col_0*
         WAMI_GMM_NUM_MODELS + k*blockDim.x + threadIdx.x] = shared_sigma[
         k*blockDim.x + threadIdx.x];
      weight[row*(WAMI_GMM_IMG_NUM_COLS*WAMI_GMM_NUM_MODELS) + col_0*
         WAMI_GMM_NUM_MODELS + k*blockDim.x + threadIdx.x] = shared_weight[
         k*blockDim.x + threadIdx.x];
   }
}

extern "C" void gpu_wami_gmm(
   u8 (foreground)[WAMI_GMM_IMG_NUM_ROWS][WAMI_GMM_IMG_NUM_COLS],
   float mu[WAMI_GMM_IMG_NUM_ROWS][WAMI_GMM_IMG_NUM_COLS][WAMI_GMM_NUM_MODELS],
   float sigma[WAMI_GMM_IMG_NUM_ROWS][WAMI_GMM_IMG_NUM_COLS]
   [WAMI_GMM_NUM_MODELS],
   float weight[WAMI_GMM_IMG_NUM_ROWS][WAMI_GMM_IMG_NUM_COLS]
   [WAMI_GMM_NUM_MODELS],
   u16 (*const frame)[WAMI_GMM_IMG_NUM_COLS])
{
   //define and allocate device arrays
   unsigned char* dev_foreground;
   float *dev_mu;
   float *dev_sigma;
   float *dev_weight;
   unsigned short *dev_frame;

   const size_t num_pixels = WAMI_GMM_IMG_NUM_ROWS * WAMI_GMM_IMG_NUM_COLS;

   CUDA_SAFE_MALLOC(dev_foreground, sizeof(unsigned char) * num_pixels);
   CUDA_SAFE_MALLOC(dev_mu, sizeof(float) * num_pixels * WAMI_GMM_NUM_MODELS);
   CUDA_SAFE_MALLOC(dev_sigma, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS);
   CUDA_SAFE_MALLOC(dev_weight, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS);
   CUDA_SAFE_MALLOC(dev_frame, sizeof(unsigned short) * num_pixels);
   CUDA_SAFE(cudaMemset(dev_foreground, 0, sizeof(unsigned char) * num_pixels));
   CUDA_SAFE(cudaMemcpy(dev_mu, mu, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS, cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_sigma, sigma, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS, cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_weight, weight, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS, cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_frame, frame, sizeof(unsigned short) * num_pixels, 
      cudaMemcpyHostToDevice));

   //compute number of blocks for kernel launch
   int num_blocks = ( WAMI_GMM_IMG_NUM_COLS + BLOCK_SIZE - 1) / BLOCK_SIZE;

   //compute number of blocks for kernel launch
   dim3 grid_dim_01(num_blocks,WAMI_GMM_IMG_NUM_ROWS,1);
   dim3 block_dim_01(BLOCK_SIZE);

   //kernel gpu_wami_gmm_kernel performs the tasks equivalent to function 
   //wami_gmm
   //launch kernel for grid_dim_01 and block_dim_01
   gpu_wami_gmm_kernel<<<grid_dim_01,block_dim_01>>>(dev_foreground, dev_mu, 
      dev_sigma, dev_weight, dev_frame);
   CUDA_SAFE(cudaDeviceSynchronize());
   //copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(foreground, dev_foreground, (sizeof(unsigned char) * 
      num_pixels), cudaMemcpyDeviceToHost));
   CUDA_SAFE(cudaMemcpy(mu, dev_mu, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS, cudaMemcpyDeviceToHost));
   CUDA_SAFE(cudaMemcpy(sigma, dev_sigma, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS, cudaMemcpyDeviceToHost));
   CUDA_SAFE(cudaMemcpy(weight, dev_weight, sizeof(float) * num_pixels * 
      WAMI_GMM_NUM_MODELS, cudaMemcpyDeviceToHost));
   //free device variables
   CUDA_SAFE(cudaFree(dev_foreground));
   CUDA_SAFE(cudaFree(dev_mu));
   CUDA_SAFE(cudaFree(dev_sigma));
   CUDA_SAFE(cudaFree(dev_weight));
   CUDA_SAFE(cudaFree(dev_frame));
}
