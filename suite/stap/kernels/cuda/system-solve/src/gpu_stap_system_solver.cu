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

#include "../lib/stap_utils.h"

/*
   define block size for kernel launch, could be variable
   tested from a value of 32 upto a maximum of 512,
   based on a different grid size, block size could be changed,
   maximum number of threads limited due to register memory
 */

#define BLOCK_SIZE 128

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { \
                     printf("CUDA CALL FAILED AT %d\n", \
                     __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE) \
                        (cudaMalloc((void**)&DP, SIZE))

__global__ void gpu_cholesky_factorization_A (complex *R)
{
   __shared__ complex shared_R[N_CHAN*TDOF][N_CHAN*TDOF];
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * 
                        blockDim.y + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N_DOP * N_BLOCKS * N_CHAN*TDOF) return;
   int dop = blockIdx.y;
   if(dop >= N_DOP) return;
   int block = blockIdx.x;
   if(block >= N_BLOCKS) return;
   int _i = threadIdx.x;
   if(_i >= N_CHAN*TDOF) return;
   int k, j, i;
   for(j = 0; j < N_CHAN*TDOF; j++)
   {
      shared_R[j][_i] = R[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*
                        (N_CHAN*TDOF*N_CHAN*TDOF)+j*(N_CHAN*TDOF)+_i];
   }
   __syncthreads();
   float Rkk_inv, Rkk_inv_sqrt;
   float2 Rkj_conj;
   float2 Rki_Rkj_conj;
   if(_i == 0)
   {
      for (k = 0; k < N_CHAN*TDOF; ++k)
      {
         /*
          * Hermitian positive definite matrices are assumed, but
          * for safety we check that the diagonal is always positive.
          */
         if( (shared_R[k][k].re <= 0))
         {
            printf("diagonal value for Hermitian matrix is <= zero for block = %d; dop = %d\n", block, dop);
            return;
         }
         /* Diagonal entries are real-valued. */
         Rkk_inv = 1.0f /shared_R[k][k].re;
         Rkk_inv_sqrt = sqrt(Rkk_inv);
         for (j = k+1; j < N_CHAN*TDOF; ++j)
         {
            Rkj_conj.x = shared_R[k][j].re;
            Rkj_conj.y = shared_R[k][j].im;
            Rkj_conj.y = -1.0f * Rkj_conj.y;
            for (i = j; i < N_CHAN*TDOF; ++i)
            { 
               Rki_Rkj_conj.x = shared_R[k][i].re * Rkj_conj.x - shared_R[k][i].im * Rkj_conj.y;
               Rki_Rkj_conj.y = shared_R[k][i].im * Rkj_conj.x + shared_R[k][i].re *  Rkj_conj.y;
               shared_R[j][i].re -= Rki_Rkj_conj.x * Rkk_inv;
               shared_R[j][i].im -= Rki_Rkj_conj.y * Rkk_inv;
            }
         }
         for (i = k; i < N_CHAN*TDOF; ++i)
         {
            shared_R[k][i].re *= Rkk_inv_sqrt;
            shared_R[k][i].im *= Rkk_inv_sqrt;
         }

      } 
   }
   complex x;
   /*
    * Copy the conjugate of the upper triangular portion of R
    * into the lower triangular portion. This is not required
    * for correctness, but can help with testing and validation
    * (e.g., correctness metrics calculated over all elements
    * will not be "diluted" by trivially correct zeros in the
    * lower diagonal region).
    */
   for (i = 0; i < N_CHAN*TDOF; ++i)
   {
      if(_i > i)
      {
         x = shared_R[i][_i];
         shared_R[_i][i].re = x.re;
         shared_R[_i][i].im = -1.0f * x.im;
      }
   }
   __syncthreads();
   for(j = 0; j < N_CHAN*TDOF; j++)
   {
      R[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*(N_CHAN*TDOF*N_CHAN*TDOF)+
        j*(N_CHAN*TDOF)+_i] = shared_R[j][_i];
   }
   __syncthreads();
}

__global__ void gpu_forward_and_back_substitution (complex *R, 
                                                   complex *steering_vectors, 
                                                   complex *x)
{
   __shared__ complex shared_x[N_CHAN*TDOF];
   __shared__ complex shared_R[N_CHAN*TDOF][N_CHAN*TDOF];

   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x + blockIdx.z*
                        gridDim.x*gridDim.y) * blockDim.x * blockDim.y
                        + (threadIdx.y * blockDim.x) + threadIdx.x;

   if(thread_number >= N_DOP * N_BLOCKS * N_STEERING * N_CHAN*TDOF) return;

   int dop = blockIdx.z;
   if(dop >= N_DOP) return;

   int block = blockIdx.y;
   if(block >= N_BLOCKS) return;

   int sv = blockIdx.x;
   if(sv >= N_STEERING) return;

   int _i = threadIdx.x;
   if(_i >= N_CHAN*TDOF) return;

   int j, k;
   float Rii_jj_inv;
   complex accum;
   complex prod;

   shared_x[_i] = x[dop*(N_BLOCKS*N_STEERING*N_CHAN*TDOF)+block*(N_STEERING*
                    N_CHAN*TDOF)+sv*(N_CHAN*TDOF)+_i];
   for(j = 0; j < N_CHAN*TDOF; j++)
   {
      shared_R[j][_i] = R[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*(N_CHAN*
                          TDOF*N_CHAN*TDOF)+j*(N_CHAN*TDOF)+_i];
   }
   __syncthreads();
   if(_i == 0 )
      for(int i = 0; i< N_CHAN*TDOF; i++)
      {
         /* First apply forward substitution */
         Rii_jj_inv = 1.0f /  shared_R[i][i].re;
         accum.re = accum.im = 0.0f;
         for (j = 0; j < i; ++j)
         {
            /*
             * Use the conjugate of the upper triangular entries
             * of R as the lower triangular entries.
             */
            prod.re = shared_R[j][i].re * shared_x[j].re + shared_R[j][i].im * 
                      shared_x[j].im;
            prod.im = shared_R[j][i].re * shared_x[j].im - shared_R[j][i].im * 
                      shared_x[j].re;

            accum.re += prod.re;
            accum.im += prod.im;
         }
         shared_x[i].re = (steering_vectors[sv*(N_CHAN*TDOF)+i].re - accum.re)*
                           Rii_jj_inv;
         shared_x[i].im = (steering_vectors[sv*(N_CHAN*TDOF)+i].im - accum.im)*
                           Rii_jj_inv;
      }
   __syncthreads();
   /* And now apply back substitution */
   if(_i == 0) 
      for (j = N_CHAN*TDOF-1; j >= 0; --j)
      {
         Rii_jj_inv = 1.0f / shared_R[j][j].re;
         accum.re = accum.im = 0.0f;
         for (k = j+1; k < N_CHAN*TDOF; ++k)
         {
            prod.re = shared_R[j][k].re * shared_x[k].re - shared_R[j][k].im * 
                      shared_x[k].im;
            prod.im = shared_R[j][k].re * shared_x[k].im + shared_R[j][k].im * 
                      shared_x[k].re;
            accum.re += prod.re;
            accum.im += prod.im;
         }
         shared_x[j].re = (shared_x[j].re - accum.re) * Rii_jj_inv;
         shared_x[j].im = (shared_x[j].im - accum.im) * Rii_jj_inv;
      }
   __syncthreads();
   x[dop*(N_BLOCKS*N_STEERING*N_CHAN*TDOF)+block*(N_STEERING*N_CHAN*TDOF)+
     sv*(N_CHAN*TDOF)+_i] = shared_x[_i];
}

extern "C" int gpu_stap_system_solver(
        complex adaptive_weights[N_DOP][N_BLOCKS][N_STEERING][N_CHAN*TDOF],
        complex (* const covariance)[N_BLOCKS][N_CHAN*TDOF][N_CHAN*TDOF],
        complex (* const steering_vectors)[N_CHAN*TDOF],
        complex cholesky_factors[N_DOP][N_BLOCKS][N_CHAN*TDOF][N_CHAN*TDOF] )
{
   complex *covariance_01 = &covariance[0][0][0][0];
   complex *steering_vectors_01 = &steering_vectors[0][0];
   const size_t num_adaptive_weight_elements = N_DOP * N_BLOCKS *
                                               N_STEERING * (N_CHAN*TDOF);
   const size_t num_covariance_elements = (TDOF*N_CHAN) * (TDOF*N_CHAN) *
                                           N_DOP * N_BLOCKS;
   const size_t num_steering_vector_elements = N_STEERING *
                                               (N_CHAN*TDOF);
   //define and allocate device arrays
   complex * dev_adaptive_weights;
   complex * dev_steering_vectors;
   complex * dev_cholesky_factors;

   CUDA_SAFE(cudaMalloc((void **)&dev_adaptive_weights, 
             num_adaptive_weight_elements * sizeof(complex)));
   CUDA_SAFE(cudaMalloc((void **)&dev_steering_vectors, 
             num_steering_vector_elements * sizeof(complex)));
   CUDA_SAFE(cudaMalloc((void **)&dev_cholesky_factors, 
             num_covariance_elements * sizeof(complex)));

   //copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemset(dev_adaptive_weights, 0, 
             num_adaptive_weight_elements * sizeof(complex)));
   CUDA_SAFE(cudaMemcpy(dev_cholesky_factors, covariance_01, 
             num_covariance_elements * sizeof(complex), 
             cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemcpy(dev_steering_vectors, steering_vectors_01, 
             num_steering_vector_elements * sizeof(complex), 
             cudaMemcpyHostToDevice));
   //compute number of blocks for kernel launch
   dim3 grid_dim(N_BLOCKS,N_DOP);
   int block_dim = N_CHAN*TDOF;
   //launch kernel for num_blocks and BLOCK_SIZE
   gpu_cholesky_factorization_A<<<grid_dim,block_dim>>>(dev_cholesky_factors);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   dim3 grid_dim2(N_STEERING,N_BLOCKS,N_DOP);
   gpu_forward_and_back_substitution<<<grid_dim2,block_dim>>>(
       dev_cholesky_factors, dev_steering_vectors, dev_adaptive_weights);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(adaptive_weights, dev_adaptive_weights, 
            (num_adaptive_weight_elements * sizeof(complex)), 
             cudaMemcpyDeviceToHost));
   //free device variables
   CUDA_SAFE(cudaFree(dev_adaptive_weights));
   CUDA_SAFE(cudaFree(dev_steering_vectors));
   CUDA_SAFE(cudaFree(dev_cholesky_factors));
   //cudaDeviceReset();
   return 0;
}
