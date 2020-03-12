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

#define PAD WAMI_DEBAYER_PAD

#define BLOCK_SIZE 128

#define CUDA_SAFE(x) if ( cudaSuccess != (x) ) { printf("CUDA CALL FAILED AT %d\n", __LINE__ ); exit(1);}
#define CUDA_SAFE_MALLOC(DP, SIZE)  (cudaMalloc((void**)&DP, SIZE))


#define PIXEL_MAX 65535

__global__ void gpu_copy_pixels_R(rgb_pixel debayered[WAMI_DEBAYER_IMG_NUM_ROWS-
   2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS]
   [WAMI_DEBAYER_IMG_NUM_COLS], int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row > (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;
   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col > (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;
   //Copy red pixels through directly
   debayered[row-PAD][col-PAD].r = bayer[row][col];
}

__global__ void gpu_copy_pixels_G( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row > (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col > (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   //Copy green pixels through directly
   debayered[row-PAD][col-PAD].g = bayer[row][col];
}

__global__ void gpu_copy_pixels_B( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y; 
   int row = row_start + row_ * 2;
   if(row > (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col > (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   //Copy blue pixels through directly
   debayered[row-PAD][col-PAD].b = bayer[row][col];
}


__global__ void gpu_interpolate_pixels_G( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;
   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   u16 pos =
      2 * bayer[row-1][col] +
      2 * bayer[row][col-1] +
      4 * bayer[row][col] +
      2 * bayer[row][col+1] +
      2 * bayer[row+1][col];
   u16 neg =
      bayer[row][col+2] +
      bayer[row-2][col] +
      bayer[row][col-2] +
      bayer[row+2][col];
   u16 pixel;
   if (pos < neg)
   {
      pixel = 0;
   }
   else
   {
      pixel = (pos - neg) >> 3;
      if (pixel > PIXEL_MAX) pixel = PIXEL_MAX; 
   }
   debayered[row-PAD][col-PAD].g = pixel;
}

__global__ void gpu_interpolate_pixels_R_at_GRB( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
int row_ = blockIdx.y;
int row = row_start + row_ * 2;
if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;
//thread ID for given bock size and number of blocks
int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
int col = col_start + col_ * 2;
if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;
const u16 pos =
((bayer[row-2][col] + bayer[row+2][col]) >> 1) +
4 * bayer[row][col-1] +
5 * bayer[row][col] +
4 * bayer[row][col+1];
const u16 neg =
bayer[row-1][col-1] +
bayer[row-1][col+1] +
bayer[row][col-2] +
bayer[row][col+2] +
bayer[row+1][col-1] +
bayer[row+1][col+1];

u16 pixel;

if (pos < neg)
{
pixel = 0;
}
else
{
pixel = (pos - neg) >> 3;
if (pixel > PIXEL_MAX) pixel = PIXEL_MAX; 
}
debayered[row-PAD][col-PAD].r = pixel;
}


__global__ void gpu_interpolate_pixels_B_at_GBR( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;
   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;
   const u16 pos =
      ((bayer[row-2][col] + bayer[row+2][col]) >> 1) +
      4 * bayer[row][col-1] +
      5 * bayer[row][col] +
      4 * bayer[row][col+1];
   const u16 neg =
      bayer[row-1][col-1] +
      bayer[row-1][col+1] +
      bayer[row][col-2] +
      bayer[row][col+2] +
      bayer[row+1][col-1] +
      bayer[row+1][col+1];
   u16 pixel;
   if (pos < neg)
   {
      pixel = 0;
   }
   else
   {
      pixel = (pos - neg) >> 3;
      if (pixel > PIXEL_MAX) pixel = PIXEL_MAX; 
   }
   debayered[row-PAD][col-PAD].b = pixel;
}

__global__ void gpu_interpolate_pixels_B( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
int row_ = blockIdx.y;
int row = row_start + row_ * 2;
if(row > (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

//thread ID for given bock size and number of blocks
int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
int col = col_start + col_ * 2;
if(col > (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

const u16 pos =
4 * bayer[row-1][col] +
((bayer[row][col-2] + bayer[row][col+2]) >> 1) +
5 * bayer[row][col] +
4 * bayer[row+1][col];
const u16 neg =
bayer[row-2][col] +
bayer[row-1][col-1] +
bayer[row-1][col+1] +
bayer[row+1][col-1] +
bayer[row+1][col+1] +
bayer[row+2][col];

u16 pixel;

if (pos < neg)
{
pixel = 0;
}
else
{
pixel = (pos - neg) >> 3;
if (pixel > PIXEL_MAX) { pixel = PIXEL_MAX; }
}
debayered[row-PAD][col-PAD].b = pixel;
}


__global__ void gpu_interpolate_pixels_R_at_GBR( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   const u16 pos =
      4 * bayer[row-1][col] +
      ((bayer[row][col-2] + bayer[row][col+2]) >> 1) +
      5 * bayer[row][col] +
      4 * bayer[row+1][col];
   const u16 neg =
      bayer[row-2][col] +
      bayer[row-1][col-1] +
      bayer[row-1][col+1] +
      bayer[row+1][col-1] +
      bayer[row+1][col+1] +
      bayer[row+2][col];

   u16 pixel;

   if (pos < neg)
   {
      pixel = 0;
   }
   else
   {
      pixel = (pos - neg) >> 3;
      if (pixel > PIXEL_MAX) pixel = PIXEL_MAX; 
   }
   debayered[row-PAD][col-PAD].r = pixel;
}

__global__ void gpu_interpolate_pixels_B_at_GRB( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   const u16 pos =
       4 * bayer[row-1][col] +
       ((bayer[row][col-2] + bayer[row][col+2]) >> 1) +
       5 * bayer[row][col] +
       4 * bayer[row+1][col];
   const u16 neg =
       bayer[row-2][col] +
       bayer[row-1][col-1] +
       bayer[row-1][col+1] +
       bayer[row+1][col-1] +
       bayer[row+1][col+1] +
       bayer[row+2][col];
   u16 pixel;
   if (pos < neg)
   {
      pixel = 0;
   }
   else
   {
      pixel = (pos - neg) >> 3;
      if (pixel > PIXEL_MAX) pixel = PIXEL_MAX; 
   }

   debayered[row-PAD][col-PAD].b = pixel;
}

__global__ void gpu_interpolate_pixels_R_at_BBB( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   u16 pos =
      2 * bayer[row-1][col-1] +
      2 * bayer[row-1][col+1] +
      6 * bayer[row][col] +
      2 * bayer[row+1][col-1] +
      2 * bayer[row+1][col+1];
   u16 neg =
      (3 * bayer[row-2][col] +
      3 * bayer[row][col-2] +
      3 * bayer[row][col+2] +
      3 * bayer[row+2][col]);
   u16 pixel;

   const u32 pos_u32 = ((u32) pos) << 1;
   const u32 neg_u32 = (u32) neg;
   if (pos_u32 < neg_u32)
   {
      pixel = 0;
   }
   else
   {
      pixel = (u16) ((pos_u32 - neg_u32) >> 4);
      if (pixel > PIXEL_MAX) pixel = PIXEL_MAX;
   }
   debayered[row-PAD][col-PAD].r = pixel;
}

__global__ void gpu_interpolate_pixels_B_at_RRR( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD], 
   u16 bayer[WAMI_DEBAYER_IMG_NUM_ROWS][WAMI_DEBAYER_IMG_NUM_COLS], 
   int row_start, int col_start)
{
   int row_ = blockIdx.y;
   int row = row_start + row_ * 2;
   if(row >= (WAMI_DEBAYER_IMG_NUM_ROWS - PAD) ) return;

   //thread ID for given bock size and number of blocks
   int col_ = (blockIdx.x * blockDim.x) + threadIdx.x;
   int col = col_start + col_ * 2;
   if(col >= (WAMI_DEBAYER_IMG_NUM_COLS - PAD) ) return;

   u16 pos =
      2 * bayer[row-1][col-1] +
      2 * bayer[row-1][col+1] +
      6 * bayer[row][col] +
      2 * bayer[row+1][col-1] +
      2 * bayer[row+1][col+1];
   u16 neg =
      (3 * bayer[row-2][col] +
      3 * bayer[row][col-2] +
      3 * bayer[row][col+2] +
      3 * bayer[row+2][col]);

   u16 pixel;

   const u32 pos_u32 = ((u32) pos) << 1;
   const u32 neg_u32 = (u32) neg;
   if (pos_u32 < neg_u32)
   {
      pixel = 0;
   }
   else
   {
      pixel = (u16) ((pos_u32 - neg_u32) >> 4);
      if (pixel > PIXEL_MAX) pixel = PIXEL_MAX;
   }
   debayered[row-PAD][col-PAD].b = pixel;
}

extern "C" void gpu_wami_debayer( rgb_pixel debayered
   [WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD][WAMI_DEBAYER_IMG_NUM_COLS-2*PAD],
   u16 (* const bayer)[WAMI_DEBAYER_IMG_NUM_COLS])
{
   const size_t num_bayer_pixels = WAMI_DEBAYER_IMG_NUM_ROWS * 
      WAMI_DEBAYER_IMG_NUM_COLS;
   const size_t num_debayer_pixels = (WAMI_DEBAYER_IMG_NUM_ROWS-2*PAD) * 
      (WAMI_DEBAYER_IMG_NUM_COLS-2*PAD);
   //define and allocate device arrays
   u16 (*dev_bayer)[WAMI_DEBAYER_IMG_NUM_COLS];
   rgb_pixel (*dev_debayer)[WAMI_DEBAYER_IMG_NUM_COLS-2*PAD];
   CUDA_SAFE_MALLOC(dev_bayer, sizeof(u16) * num_bayer_pixels);
   CUDA_SAFE_MALLOC(dev_debayer, sizeof(rgb_pixel) * num_debayer_pixels);
   CUDA_SAFE(cudaMemcpy(dev_bayer, bayer, sizeof(u16) * num_bayer_pixels, 
      cudaMemcpyHostToDevice));
   CUDA_SAFE(cudaMemset(dev_debayer, 0, sizeof(rgb_pixel) * 
      num_debayer_pixels));
   //compute number of blocks for kernel launch
   int num_blocks_x = ( (WAMI_DEBAYER_IMG_NUM_COLS - 2*PAD)/2 + BLOCK_SIZE - 1)
      / BLOCK_SIZE;
   int num_blocks_y =  (WAMI_DEBAYER_IMG_NUM_ROWS - 2*PAD)/2;
   //compute number of blocks for kernel launch
   dim3 grid_dim_01(num_blocks_x,num_blocks_y,1);
   dim3 block_dim_01(BLOCK_SIZE);
   //Copy red pixels through directly
   gpu_copy_pixels_R<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD, PAD);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Copy top-right green pixels through directly
   gpu_copy_pixels_G<<<grid_dim_01,block_dim_01>>>(dev_debayer, dev_bayer, 
      PAD, PAD+1);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Copy bottom-left green pixels through directly
   gpu_copy_pixels_G<<<grid_dim_01,block_dim_01>>>(dev_debayer, dev_bayer, 
      PAD+1, PAD);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Copy blue pixels through directly
   gpu_copy_pixels_B<<<grid_dim_01,block_dim_01>>>(dev_debayer, dev_bayer, 
      PAD+1, PAD+1);
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate green pixels at red pixels
   gpu_interpolate_pixels_G<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD, PAD );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate green pixels at blue pixels
   gpu_interpolate_pixels_G<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD+1, PAD+1 );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate red pixels at green pixels, red row, blue column
   gpu_interpolate_pixels_R_at_GRB<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD, PAD+1 );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate blue pixels at green pixels, blue row, red column
   gpu_interpolate_pixels_B_at_GBR<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD+1, PAD );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate red pixels at green pixels, blue row, red column
   gpu_interpolate_pixels_R_at_GBR<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD+1, PAD );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate blue pixels at green pixels, red row, blue column
   gpu_interpolate_pixels_B_at_GRB<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD, PAD+1 );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate red pixels at blue pixels, blue row, blue column
   gpu_interpolate_pixels_R_at_BBB<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD+1, PAD+1 );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //Interpolate blue pixels at red pixels, red row, red column
   gpu_interpolate_pixels_B_at_RRR<<<grid_dim_01,block_dim_01>>>(dev_debayer, 
      dev_bayer, PAD, PAD );
   //synchronize threads
   CUDA_SAFE(cudaDeviceSynchronize());
   //copy variable out back from Device to Host
   CUDA_SAFE(cudaMemcpy(debayered, dev_debayer, sizeof(rgb_pixel) * 
      num_debayer_pixels, cudaMemcpyDeviceToHost));
   //free device variables
   cudaFree(dev_bayer);
   cudaFree(dev_debayer);
}
