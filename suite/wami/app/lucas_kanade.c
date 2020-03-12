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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "wami_lucas_kanade.h"
#include "timer.h"
#include "matrix_ops.h"
#include "gradient.h"

/* macro for (ideally) clean exits on failure */
#define K2_FAIL_FREE_EXIT() \
    if(error_img)  free (error_img); \
    if(img_dx)     free (img_dx); \
    if(img_dy)     free (img_dy); \
    if(nabla_Ix)   free (nabla_Ix); \
    if(nabla_Iy)   free (nabla_Iy); \
    if(I_steepest) free (I_steepest); \
    if(sd_delta_p) free (sd_delta_p); \
    if(delta_p)    free (delta_p); \
    if(H)          free (H); \
    if(H_wsp)      free (H_wsp); \
    if(H_inv)      free (H_inv); \
    return -1;

/* Registers tmplt to img and produces a warped output in warped_tmplt to align 
 * the images.  Note that the warp parameters are not zeroed such that you
 * can sart with the previous settings or reinitialize at each step.
 * RETURNS non-zero on error
 */
int
kernel2_lucas_kanade (
  fltPixel_t * img, 
  fltPixel_t * tmplt, 
  float      * affine_warp, 
  int num_iterations,
  int nRows, int nCols,
  fltPixel_t * warped_tmplt)
{
  /* Variables --------------------------------------------------------------- */

  /* Input images */
  fltPixel_t * error_img  = NULL;
  /* Gradients of input image */
  fltPixel_t * img_dx     = NULL;
  fltPixel_t * img_dy     = NULL;
  /* Warped input gradients */
  fltPixel_t * nabla_Ix   = NULL;
  fltPixel_t * nabla_Iy   = NULL;
  /* Steepest descent */
  fltPixel_t * I_steepest = NULL;
  fltPixel_t * sd_delta_p = NULL;
  fltPixel_t * delta_p    = NULL;
  /* Hessian */
  float * H               = NULL;
  float * H_wsp           = NULL;
  float * H_inv           = NULL;

  /* Alias for consistency with original matlab code */
  float * IWxp = warped_tmplt;

  int iteration = 0;
  int N_p = 6;

  PRINT_STAT_STRING ("kernel", "lucas_kanade");
  PRINT_STAT_INT ("rows", nRows);
  PRINT_STAT_INT ("columns", nCols);

  /* Allocate and check ------------------------------------------------------ */
  error_img  = calloc (nRows * nCols, sizeof(fltPixel_t));
  img_dx     = calloc (nRows * nCols, sizeof(fltPixel_t));
  img_dy     = calloc (nRows * nCols, sizeof(fltPixel_t));
  nabla_Ix   = calloc (nRows * nCols, sizeof(fltPixel_t));
  nabla_Iy   = calloc (nRows * nCols, sizeof(fltPixel_t));
  I_steepest = calloc (6 * nRows * nCols, sizeof(fltPixel_t));
  sd_delta_p = calloc (6 * 1, sizeof(fltPixel_t));
  delta_p    = calloc (6 * 1, sizeof(fltPixel_t));
  H          = calloc (6 * 6, sizeof(float));
  H_wsp      = calloc (6 * 6, sizeof(float));
  H_inv      = calloc (6 * 6, sizeof(float));

  if (!error_img || !img_dx || !img_dy || !nabla_Ix || !nabla_Iy || 
      !I_steepest || !sd_delta_p || !delta_p || !H || !H_wsp || !H_inv) {

    fprintf(stderr, "ERROR: Allocation failed.\n");
    K2_FAIL_FREE_EXIT();
  }

  /* Start computation -------------------------------------------------------- */
  /* From here, this will follow the Matlab forward-additive implementation from 
    Lucas-Kanade 20 Years On
    Iain Matthews, Simon Baker, Carnegie Mellon University, Pittsburgh
    $Id$
  */
  
  /* Pre-computable things --------------------------------------------------- */

  /* 3a) Compute gradients */
  if(gradientXY(img, nCols, nRows, img_dx, img_dy)) {
    fprintf(stderr, "ERROR: gradient failed.\n");
    K2_FAIL_FREE_EXIT();
  }

  /* 4) Jacobian (defined in steepest descent step) */
  /* nothing to do...                               */

  /* p - affine parameters */
  /* warp_p[0] =  horizontal compression */
  /* warp_p[1] =  horizontal distortion */
  /* warp_p[2] =  vertical distortion */
  /* warp_p[3] =  vertical compression */
  /* warp_p[4] =  horizontal translation */
  /* warp_p[5] =  vertical translation */

  /* Main loop --------------------------------------------------------------- */
  for(iteration = 0; iteration < num_iterations; iteration++) {
    printf(",\n\t\"warp matrix\": [[%8e, %8e, %8e], [%8e, %8e, %8e]]", 
      affine_warp[0], affine_warp[1], affine_warp[2], 
      affine_warp[3], affine_warp[4], affine_warp[5]);

    /* 1) Warp image with current parameters */
    tic ();
    warp_image(img, nRows, nCols, affine_warp, IWxp);
    PRINT_STAT_COUNT_DOUBLE (iteration, "1)time_warp", toc ());

    /* 2) compute error image */
    tic ();
    if(subtract(tmplt, IWxp, nRows, nCols, error_img)) {
      fprintf(stderr, "ERROR: Subtraction failed.\n");
      K2_FAIL_FREE_EXIT();
    }
    PRINT_STAT_COUNT_DOUBLE (iteration, "2)time_subtract", toc ());

    /* 3b) Warp the gradient I with W(x;p) */
    tic ();
    warp_image (img_dx, nRows, nCols, affine_warp, nabla_Ix);
    warp_image (img_dy, nRows, nCols, affine_warp, nabla_Iy);
    PRINT_STAT_COUNT_DOUBLE (iteration, "3b)time_warp", toc ());

    /* 4) Jacobian and ...
       5) Compute the steepest descent images Gradient * Jacobian */
    tic ();
    steepest_descent (nabla_Ix, nabla_Iy, nRows, nCols, I_steepest);
    PRINT_STAT_COUNT_DOUBLE (iteration, "4+5)time_steepest_descent", toc ());

    /* 6) Compute the Hessian matrix and inverse */
    tic ();
    hessian (I_steepest, nRows, nCols, N_p, H);
    invert_gauss_jordan(H, H_wsp, 6, 6, H_inv);
    PRINT_STAT_COUNT_DOUBLE (iteration, "6)time_hessian", toc ());

    /* 7) Update based on steepest descent */
    tic ();
    sd_update(I_steepest, error_img, N_p, nRows, nCols, sd_delta_p);
    PRINT_STAT_COUNT_DOUBLE (iteration, "7)time_sd_update", toc ());

    /* 8) compute gradient descent parameter updates */
    tic ();
    mult(H_inv, sd_delta_p, 6, 1, 6, delta_p);
    PRINT_STAT_COUNT_DOUBLE (iteration, "8)time_sd_param_updates", toc ());

    /* reuse sd_delta_p to store the reshaped delta_p */
    reshape(delta_p, 6, 1, 2, 3, sd_delta_p);

    /* 9) update warp parameters */
    tic ();
    add(affine_warp, sd_delta_p, 2, 3, affine_warp);
    PRINT_STAT_COUNT_DOUBLE (iteration, "9)time_warp_param_updates", toc ());
  }

  tic ();
  warp_image(img, nRows, nCols, affine_warp, IWxp);
  PRINT_STAT_DOUBLE ("final)time_warp", toc ());

  if(error_img)  free (error_img); 
  if(img_dx)     free (img_dx); 
  if(img_dy)     free (img_dy); 
  if(nabla_Ix)   free (nabla_Ix); 
  if(nabla_Iy)   free (nabla_Iy); 
  if(I_steepest) free (I_steepest); 
  if(sd_delta_p) free (sd_delta_p); 
  if(delta_p)    free (delta_p); 
  if(H)          free (H); 
  if(H_wsp)      free (H_wsp); 
  if(H_inv)      free (H_inv); 

  return 0;
}

int
sd_update(fltPixel_t * VI_dW_dp, fltPixel_t * error_img, int N_p, int nRows, int nCols, fltPixel_t * sd_delta_p)
{
  int i, x, y;

  if(!VI_dW_dp || !error_img || !sd_delta_p)
    return -1;

  if(N_p < 0 || nRows < 0 || nCols < 0)
    return -2;

  if(N_p == 0 || nRows == 0 || nCols == 0)
    return 0;

  for(i = 0; i < N_p; i++) {
    sd_delta_p[i] = 0.0;
  }

  for(i = 0; i < N_p; i++) {
    for (y = 0; y < nRows; y++) {
      for (x = 0; x < nCols; x++) {
	sd_delta_p[i] += error_img[y * nCols + x] * VI_dW_dp[y * N_p * nCols + (x + i * nCols)];
      }
    }
  }

  return 0;
}

void
warp_image (fltPixel_t *Iin, int nCols, int nRows, float *W_xp, fltPixel_t *Iout)
{
  int x, y;
  float Tlocalx, Tlocaly;
  float compa0, compa1, compb0, compb1;
  int index = 0;

  compb0 = W_xp[2];
  compb1 = W_xp[5];

  for (y = 0; y < nRows; y++) {
    compa0 = W_xp[1] * ((float) y) + compb0;
    compa1 = (1.0 + W_xp[4]) * ((float) y) + compb1;

    for (x = 0; x < nCols; x++) {
      Tlocalx = (1.0 + W_xp[0]) * ((float) x) + compa0;
      Tlocaly = W_xp[3] * ((float) x) + compa1;

      Iout[index] = interpolate (Tlocalx, Tlocaly, nCols, nRows, Iin);
      index++;
    }
  }

}

void
steepest_descent (fltPixel_t *gradX_warped, fltPixel_t *gradY_warped, int nCols, int nRows, fltPixel_t *I_steepest)
{
  int k;
  int x, y;
  float Jacobian_x[6], Jacobian_y[6];
  int index, j_index;

  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      index = y * nCols + x;

      Jacobian_x[0] = (float) x;
      Jacobian_x[1] = 0.0;
      Jacobian_x[2] = (float) y;
      Jacobian_x[3] = 0.0;
      Jacobian_x[4] = 1.0;
      Jacobian_x[5] = 0.0;

      Jacobian_y[0] = 0.0;
      Jacobian_y[1] = (float) x;
      Jacobian_y[2] = 0.0;
      Jacobian_y[3] = (float) y;
      Jacobian_y[4] = 0.0;
      Jacobian_y[5] = 1.0;

      for (k = 0; k < 6; k++) {
	j_index = (6 * y * nCols) + (nCols * k) + x;
	I_steepest[j_index] = (Jacobian_x[k] * gradX_warped[index]) + (Jacobian_y[k] * gradY_warped[index]);
      }
    }
  }
}

void
hessian (fltPixel_t *I_steepest, int nCols, int nRows, int np, float *H)
{
  int i, j;
  int x, y;

  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      H[i * np + j] = 0;
    }
  }

  /* compare each image in the 6-wide I_steepest to each other image */
  for (i = 0; i < np; i++) {
    for (j = 0; j < np; j++) {
      
      /* sum the element-wise product of the images */
      double total = 0.0;
      for (y = 0; y < nRows; y++) {
	for (x = 0; x < nCols; x++) {
	  int index1 = (np * y * nCols) + (nCols * i) + x;
	  int index2 = (np * y * nCols) + (nCols * j) + x;
	  total += ((double)I_steepest[index1]) * ((double)I_steepest[index2]);
	}
      }

      H[np * j + i] = total;
    }
  }

}

