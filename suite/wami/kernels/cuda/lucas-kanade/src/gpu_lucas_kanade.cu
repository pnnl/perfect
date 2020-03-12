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
#include "wami_lucas_kanade.h"

#define BLOCK_SIZE 32

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

__global__ void gpu_warp_image (float *Iin, int nCols, int nRows, float *W_xp_,    float *Iout)
{
   __shared__ float W_xp[6];

   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * 
      blockDim.y + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= nRows * nCols) return;
   int y_ = blockIdx.y;
   int y = y_ + 1;
   if(y > nRows) return;
   //thread ID for given bock size and number of blocks
   int x_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int x = x_ + 1;
   if(x > nCols) return;

   if(threadIdx.x == 0) 
   {
      W_xp[0] = W_xp_[0];
      W_xp[1] = W_xp_[1];
      W_xp[2] = W_xp_[2];
      W_xp[3] = W_xp_[3];
      W_xp[4] = W_xp_[4];
      W_xp[5] = W_xp_[5];
   } 
   __syncthreads();
   float compb0, compb1;
   compb0 = W_xp[2];
   compb1 = W_xp[5];
   float compa0, compa1;
   compa0 = W_xp[1] * ((float) y) + compb0;
   compa1 = W_xp[4] * ((float) y) + compb1;
   float Tlocalx, Tlocaly;
   Tlocalx = W_xp[0] * ((float) x) + compa0;
   Tlocaly = W_xp[3] * ((float) x) + compa1;
   float interpolate_val;
   // Linear interpolation variables 
   int xBas0, xBas1, yBas0, yBas1;
   float perc[4] = {0, 0, 0, 0};
   float xCom, yCom, xComi, yComi;
   float color[4] = {0, 0, 0, 0};
   // Rounded location 
   float fTlocalx, fTlocaly;
   // Determine the coordinates of the pixel(s) which will become the 
   // current pixel (using linear interpolation)
   fTlocalx = floor(Tlocalx);
   fTlocaly = floor(Tlocaly);
   xBas0 = (int) fTlocalx;
   yBas0 = (int) fTlocaly;
   xBas1 = xBas0 + 1;
   yBas1 = yBas0 + 1;
   // Linear interpolation constants (percentages)
   xCom = Tlocalx - fTlocalx;
   yCom = Tlocaly - fTlocaly;
   xComi = (1.0f - xCom);
   yComi = (1.0f - yCom);
   perc[0] = xComi * yComi;
   perc[1] = xComi * yCom;
   perc[2] = xCom * yComi;
   perc[3] = xCom * yCom;

   if (xBas0 < 0) {
      Iout[(y-1)*nCols+ (x-1)] = 0.0f;
      return;
   }

   if (yBas0 < 0) {
      Iout[(y-1)*nCols+ (x-1)] = 0.0f;
      return;
   }

   if (xBas1 > (nCols - 1)) {
      Iout[(y-1)*nCols+ (x-1)] = 0.0f;
      return;
   }

   if (yBas1 > (nRows - 1)) {
      Iout[(y-1)*nCols+ (x-1)] = 0.0f;
      return;
   }

   color[0] = Iin[yBas0 * nCols + xBas0];
   color[1] = Iin[yBas1 * nCols + xBas0];
   color[2] = Iin[yBas0 * nCols + xBas1];
   color[3] = Iin[yBas1 * nCols + xBas1];

   interpolate_val = 
      color[0] * perc[0]
      + color[1] * perc[1]
      + color[2] * perc[2]
      + color[3] * perc[3];

   Iout[(y-1)*nCols+ (x-1)] = interpolate_val;
}

__global__ void gpu_steepest_descent (float *gradX_warped, float *gradY_warped, 
   int nCols, int nRows, float *I_steepest)
{
   __shared__ float Jacobian_x[6][BLOCK_SIZE];
   __shared__ float Jacobian_y[6][BLOCK_SIZE];

   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * 
      blockDim.y + (threadIdx.y * blockDim.x) + threadIdx.x;

   if(thread_number >= nRows * nCols) return;

   //thread ID for given bock size and number of blocks
   int x = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(x >= nCols) return;

   int y = blockIdx.y;
   if(y >= nRows) return;

   Jacobian_x[0][threadIdx.x] = (float) x;
   Jacobian_x[1][threadIdx.x] = 0.0f;
   Jacobian_x[2][threadIdx.x] = (float) y;
   Jacobian_x[3][threadIdx.x] = 0.0f;
   Jacobian_x[4][threadIdx.x] = 1.0f;
   Jacobian_x[5][threadIdx.x] = 0.0f;

   Jacobian_y[0][threadIdx.x] = 0.0f;
   Jacobian_y[1][threadIdx.x] = (float) x;
   Jacobian_y[2][threadIdx.x] = 0.0f;
   Jacobian_y[3][threadIdx.x] = (float) y;
   Jacobian_y[4][threadIdx.x] = 0.0f;
   Jacobian_y[5][threadIdx.x] = 1.0f;

   int index, j_index;
   int k;
   index = y * nCols + x;

   for (k = 0; k < 6; k++) 
   {
      j_index = (6 * y * nCols) + (nCols * k) + x;
      I_steepest[j_index] = (Jacobian_x[k][threadIdx.x] * gradX_warped[index]) 
         + (Jacobian_y[k][threadIdx.x] * gradY_warped[index]);
   }
}

/* Atomic double operation */
__device__ double double_add(double *address, double val){
   unsigned long long int *address_l = (unsigned long long *) address;
   unsigned long long old = *address_l , assumed;
   do{
      assumed = old;
      old = atomicCAS(address_l, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
   } while (assumed != old);
   return __longlong_as_double(old);
}

__global__ void gpu_hessian (float *I_steepest, int nCols, int nRows, int np,
    double *H)
{
   __shared__ double shared_H[6*6][BLOCK_SIZE];
   int i,j;
   for (i = 0; i < np; i++) {
      for (j = 0; j < np; j++) {
         shared_H[6*i + j][threadIdx.x] = 0.0f;
      }
   }
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * 
      blockDim.y + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= nRows * nCols) return;
   //thread ID for given bock size and number of blocks
   int x = (blockIdx.x * blockDim.x) + threadIdx.x;
   if(x >= nCols) return;

   int y = blockIdx.y;
   if(y >= nRows) return;

   for (i = 0; i < np; i++) {
      for (j = 0; j < np; j++) {
         int index1 = (6 * y * nCols) + (nCols * i) + x;
         int index2 = (6 * y * nCols) + (nCols * j) + x;
         shared_H[6*i + j][threadIdx.x] = I_steepest[index1] * 
            I_steepest[index2];
      }
   }
   __syncthreads();

   //Perform reduction on shared_H
   for(int offset = blockDim.x>>1; offset > 0; offset >>=1)
   {
      if(threadIdx.x < offset)
      {
         for (i = 0; i < np; i++) {
            for (j = 0; j < np; j++) {
               shared_H[6*i + j][threadIdx.x] +=shared_H[6*i + j]
                  [threadIdx.x + offset];
               shared_H[6*i + j][threadIdx.x + offset] = 0.0f;
            }
         }
      }
      __syncthreads();
   }

/* Numerically Unstable */
   if(threadIdx.x == 0)
   {

      for (i = 0; i < np; i++) {
         for (j = 0; j < np; j++) {
            //atomicAdd(&H[6*i + j], shared_H[6*i + j][threadIdx.x]);
            double_add(&H[6*i + j], shared_H[6*i + j][threadIdx.x]);
         }
      }
   }
}

extern "C" void gpu_lucas_kanade(fltPixel_t *gradX, fltPixel_t *gradY,
   int nCols, int nRows, double *H) 
{
   int N = nCols;
   int M = nRows;
   /* Parameter set p */
   float p_in[6];
   /* Warp/ pixel transformation matrix */
   float W_xp[9];
   /* p */
   p_in[0] = -0.035f;      /* horizontal compression */
   p_in[1] = 0.01f;        /* horizontal distortion */
   p_in[2] = -0.01f;       /* vertical distortion */
   p_in[3] = -0.035f;      /* vertical compression */
   p_in[4] = 5.5f;         /* horizontal translation */
   p_in[5] = 5.5f;         /* vertical translation */

   /* W(x;p) */
   W_xp[0] = 1.0f + p_in[0];  W_xp[1] = p_in[2];        W_xp[2] = p_in[4];
   W_xp[3] = p_in[1];        W_xp[4] = 1.0f + p_in[3];  W_xp[5] = p_in[5];
   W_xp[6] = 0.0f;            W_xp[7] = 0.0f;            W_xp[8] = 1.0f;

   //define and allocate device arrays
   float *dev_W_xp;

   fltPixel_t *dev_gradX; 
   fltPixel_t *dev_gradY; 

   fltPixel_t *dev_gradX_warped; 
   fltPixel_t *dev_gradY_warped; 

   fltPixel_t *dev_I_steepest;
   double  *dev_H;

   CUDA_SAFE_MALLOC(dev_W_xp, sizeof(float) * 9);
   CUDA_SAFE_MALLOC(dev_gradX, sizeof(fltPixel_t) * M * N);
   CUDA_SAFE_MALLOC(dev_gradY, sizeof(fltPixel_t) * M * N);

   CUDA_SAFE_MALLOC(dev_gradX_warped, sizeof(fltPixel_t) * M * N);
   CUDA_SAFE_MALLOC(dev_gradY_warped, sizeof(fltPixel_t) * M * N);

   CUDA_SAFE_MALLOC(dev_I_steepest, sizeof(fltPixel_t) * 6 *  M * N);

   CUDA_SAFE_MALLOC(dev_H, sizeof(double) * 6 *  6);

   CUDA_SAFE(cudaMemcpy(dev_W_xp, W_xp, sizeof(float) * 9, 
      cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_gradX, gradX, sizeof(fltPixel_t) * M * N, 
      cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_gradY, gradY, sizeof(fltPixel_t) * M * N, 
      cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemset(dev_gradX_warped, 0, sizeof(fltPixel_t) * M * N));
   CUDA_SAFE(cudaMemset(dev_gradY_warped, 0, sizeof(fltPixel_t) * M * N));
   CUDA_SAFE(cudaMemset(dev_I_steepest, 0, sizeof(fltPixel_t) * 6 * M * N));
   CUDA_SAFE(cudaMemset(dev_H, 0, sizeof(double) * 6 * 6));
   //compute number of blocks for kernel launch
   int num_blocks = ( nCols + BLOCK_SIZE - 1) / BLOCK_SIZE;
   //compute number of blocks for kernel launch
   dim3 grid_dim_01(num_blocks,nRows,1);
   dim3 block_dim_01(BLOCK_SIZE);
   gpu_warp_image<<<grid_dim_01,block_dim_01>>>(dev_gradX, nCols, nRows, 
      dev_W_xp, dev_gradX_warped);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   gpu_warp_image<<<grid_dim_01,block_dim_01>>>(dev_gradY, nCols, nRows, 
      dev_W_xp, dev_gradY_warped);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   gpu_steepest_descent<<<grid_dim_01,block_dim_01>>> (dev_gradX_warped, 
      dev_gradY_warped, nCols, nRows, dev_I_steepest);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   gpu_hessian<<<grid_dim_01,block_dim_01>>> (dev_I_steepest, nCols, nRows, 6, 
      dev_H);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(H, dev_H, sizeof(double) * 6 *  6, 
      cudaMemcpyDeviceToHost));
   //free device variables
   cudaFree(dev_W_xp);
   cudaFree(dev_gradX);
   cudaFree(dev_gradY);
   cudaFree(dev_gradX_warped);
   cudaFree(dev_gradY_warped);
   cudaFree(dev_I_steepest);
   cudaFree(dev_H);
}
