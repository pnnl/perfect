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
   tested from a value of 64 upto a maximum of 1024
   based on a different grid size, block size could be changed
 */

#define BLOCK_SIZE 1

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { \
    printf("CUDA CALL FAILED AT %d\n", __LINE__  \
    ); exit(1);}

#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))

__global__ void gpu_extract_snapshot_a (complex *datacube, complex *covariance)
{
   __shared__ complex shared_snapshot[N_CHAN*TDOF];
   __shared__ complex shared_conj_snapshot[N_CHAN*TDOF];
   // Thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * blockDim.y 
                        + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N_DOP * N_BLOCKS * N_CHAN * TDOF) return;
   int block = blockIdx.x;
   if(block >= N_BLOCKS) return;
   int dop_index = blockIdx.y;
   if(dop_index >= N_DOP) return;
   if(threadIdx.x >= N_CHAN * TDOF) return;
   int first_cell = block*TRAINING_BLOCK_SIZE;
   int last_cell = (block+1)*TRAINING_BLOCK_SIZE-1;
   for (int cell = first_cell; cell <= last_cell; ++cell)
   {
      int dof;
      int chan = threadIdx.x;
      if(chan < N_CHAN) 
      {
         for (dof = 0; dof < TDOF; ++dof)
         {
            int dop = dop_index - (TDOF-1)/2 + dof;
            if (dop < 0) { dop += N_DOP; }
            if (dop >= N_DOP) { dop -= N_DOP; } 
            shared_snapshot[chan*TDOF+dof] = datacube[chan*(N_DOP*N_RANGE)+
                                             dop*N_RANGE+cell];
            // calculate conjugate of snapshot
            shared_conj_snapshot[chan*TDOF+dof] = shared_snapshot[chan*TDOF+
                                                  dof];
	    shared_conj_snapshot[chan*TDOF+dof].im = -1.0f * shared_snapshot
                                                     [chan*TDOF+dof].im;
         }
      }
      __syncthreads();
      int i;
      int j = threadIdx.x;
      complex x;
      for (i = 0; i < N_CHAN*TDOF; ++i)
      {
         /* Exploit conjugate symmetry by only accumulating along
          * the diagonal and below. */
         if(j <= threadIdx.x)
         {
            x.re = shared_snapshot[i].re * shared_conj_snapshot[j].re - 
                   shared_snapshot[i].im * shared_conj_snapshot[j].im;
            x.im = shared_snapshot[i].re * shared_conj_snapshot[j].im + 
                   shared_snapshot[i].im * shared_conj_snapshot[j].re;
            atomicAdd(&covariance[dop_index*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+
                      block*(N_CHAN*TDOF*N_CHAN*TDOF)+ i*(N_CHAN*TDOF)+j].re, 
                      x.re);
            atomicAdd(&covariance[dop_index*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+
                      block*(N_CHAN*TDOF*N_CHAN*TDOF)+i*(N_CHAN*TDOF)+j].im, 
                      x.im);
         }
      }
   }
}

__global__ void gpu_copy_triangular (complex *covariance)
{
   //thread ID for given bock size and number of blocks
   int thread_number = (blockIdx.x + blockIdx.y*gridDim.x) * blockDim.x * 
                        blockDim.y + (threadIdx.y * blockDim.x) + threadIdx.x;
   if(thread_number >= N_DOP * N_BLOCKS * N_CHAN * TDOF) return;
   int block = blockIdx.x;
   if(block >= N_BLOCKS) return;
   int dop = blockIdx.y;
   if(dop >= N_DOP) return;
   int j = threadIdx.x;
   if(j  >= N_CHAN * TDOF) return;
   int i;
   complex x;
   /*
    * The covariance matrices are conjugate symmetric, so
    * we copy the conjugate of the lower triangular portion
    * into the upper triangular portion.
    */
   for (i = 0; i < N_CHAN*TDOF; ++i)
   {
      if(j > i)
      {
         x = covariance[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*
                       (N_CHAN*TDOF*N_CHAN*TDOF)+j*(N_CHAN*TDOF)+i];
         covariance[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*
                    (N_CHAN*TDOF*N_CHAN*TDOF)+i*(N_CHAN*TDOF)+j].re = x.re;
         covariance[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*
                    (N_CHAN*TDOF*N_CHAN*TDOF)+i*(N_CHAN*TDOF)+j].im = 
                    -1.0f * x.im;
      }
   }
   /*
    * Normalize the covariance matrices by dividing by the
    * number of training samples.
    */
   for (i = 0; i < N_CHAN*TDOF; ++i)
   {
      covariance[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*
                (N_CHAN*TDOF*N_CHAN*TDOF)+i*(N_CHAN*TDOF)+j].re *= 
                (1.0f/TRAINING_BLOCK_SIZE);
      covariance[dop*(N_BLOCKS*N_CHAN*TDOF*N_CHAN*TDOF)+block*
                (N_CHAN*TDOF*N_CHAN*TDOF)+i*(N_CHAN*TDOF)+j].im *= 
                (1.0f/TRAINING_BLOCK_SIZE);
   }
}

extern "C" int gpu_stap_compute_covariance_estimate(
   complex covariance[N_DOP][N_BLOCKS][N_STEERING][N_CHAN*TDOF],
   complex (* const datacube)[N_DOP][N_RANGE])
{
   /*
    * It is assumed for simplicity that the training block
    * size evenly divides the range swath.
    */
   assert(N_RANGE % TRAINING_BLOCK_SIZE == 0);
   //define and allocate device arrays
   complex *dev_datacube;
   complex *dev_covariance;
   const int num_datacube = N_CHAN * N_DOP * N_RANGE; 
   const int num_covariance_elements = N_DOP * N_BLOCKS * (TDOF*N_CHAN) * (TDOF*N_CHAN); 

   CUDA_SAFE(cudaMalloc((void **)&dev_datacube, num_datacube * sizeof(complex)));
   CUDA_SAFE(cudaMalloc((void **)&dev_covariance, num_covariance_elements * sizeof(complex)));

   //copy data from host to device or initialize device variables
   CUDA_SAFE(cudaMemcpy(dev_datacube, datacube, num_datacube * sizeof(complex), cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemset(dev_covariance, 0, num_covariance_elements * sizeof(complex)));

   /* kernels gpu_extract_snapshot_a and gpu_copy_triangular are launched N_BLOCKS number of times */

   //compute number of blocks for kernel launch
   dim3 grid_dim_01(N_BLOCKS,N_DOP,1);
   dim3 block_dim_01(N_CHAN * TDOF);

   //launch kernel for num_blocks and BLOCK_SIZE
   gpu_extract_snapshot_a<<<grid_dim_01,block_dim_01>>>(dev_datacube, dev_covariance);

   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());

   //launch kernel for num_blocks and BLOCK_SIZE
   gpu_copy_triangular<<<grid_dim_01,block_dim_01>>>(dev_covariance);

   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   
   //copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(covariance, dev_covariance, (num_covariance_elements * sizeof(complex)), cudaMemcpyDeviceToHost));

   //free device variables
   CUDA_SAFE(cudaFree(dev_datacube));
   CUDA_SAFE(cudaFree(dev_covariance));

   return 0;
}








