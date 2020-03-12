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
#include <time.h>

#include "octave_u16.h"
#include "rgb_to_grayscale.h"
#include "fileio.h"

/* from kernel 1 - debayer */
#include "wami_params.h"
#include "wami_utils.h"
#include "wami_debayer.h"

/* from kernel 2 - lucas-kanade */
#include "timer.h"
#include "octave.h"
#include "matrix_ops.h"
#include "wami_lucas_kanade.h"

/* from kernel 3 - GMM Change Detection */
#include "wami_gmm.h"

/* static configuration */
#define ITERATIONS 20

#define ENABLE_CORRECTNESS_CHECKING
#define WRITE_OUTPUT_TO_DISK

/* disagreement after change detection up to 1% of total pixels in image (MxN) */
#define ERROR_PERCENTAGE 0.01f

#if INPUT_SIZE == INPUT_SIZE_SMALL
    static const char *output_filename = "../inout/small_app_output.m";
    static const char *input_filename = "../inout/small_app_input.bin";
#elif INPUT_SIZE == INPUT_SIZE_MEDIUM
    static const char *output_filename = "../inout/medium_app_output.m";
    static const char *input_filename = "../inout/medium_app_input.bin";
#elif INPUT_SIZE == INPUT_SIZE_LARGE
    static const char *output_filename = "../inout/large_app_output.m";
    static const char *input_filename = "../inout/large_app_input.bin";
#else
    #error "Unhandled value for INPUT_SIZE"
#endif

/* local functions */
int build_test(const char * outfile, int start, int stop, int argc, char * argv[]);
int run_test(const char * infile, const char * outfile);

int
main(int argc, char * argv[]) {
  if (getenv("WAMI_BUILD_TEST")) {
    char outfile[256];

    if(argc < 6) {
      printf("Test Building Mode\n"
             "Usage: %s <filename> <start_frame> <stop_frame> <img1.mat> <img2.mat> ...\n\n"
	     "  Note that all images must have the same dimensions.\n"
	     "  Start frame and stop frame are the counts where the state, inputs, and\n"
	     "  outputs will be captured from.\n", argv[0]);
      return -1;
    }
    build_test(argv[1], atoi(argv[2]), atoi(argv[3]), argc - 3, argv + 3);
    /* C89 lacks snprintf, so length-check the strings before calling snprintf. */
    if (strlen(argv[1]) + strlen("_testing_output") + 1 > sizeof(outfile))
    {
        printf("Error: outfile filename is too large (%lu > %lu)\n",
            strlen(argv[1])+strlen("_testing_output")+1, sizeof(outfile));
        return -1;
    }
    sprintf(outfile, "%s_testing_output", argv[1]);
    run_test(argv[1], outfile);
    return 0;
  }

  if(argc > 1 && (!strncmp(argv[1], "-h", 2) || !strncmp(argv[1], "-?", 2) || !strncmp(argv[1], "help", 4))) {
    printf("Usage: build a new test: set the environment variable WAMI_BUILD_TEST\n"
           "       run a test:       pass the test file and output file on the command\n" 
	   "                         line or run without one to use the default:\n"
	   "input : %s\noutput: %s\n", input_filename, output_filename);
    return -1;
  }

  return run_test(argc > 2 ? argv[1] : input_filename, argc > 2 ? argv[2] : output_filename);
}

char *
get_name(char * name, int count) {
  static char str[256];
  const size_t INT_LEN_AND_NULL_TERMINATOR = 12;
  /* C89 lacks snprintf, so length-check the strings before calling snprintf. */
  if (strlen(name) > sizeof(str) - INT_LEN_AND_NULL_TERMINATOR)
  {
    printf("%s:%d: Name prefix in get_name too long.\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }
  sprintf(str, "%s%d", name, count);
  return str;
}

int
run_test(const char * infile, const char * outfile)
{
  FILE * outfp = NULL;

  /* Input images */
  u32          nTestImages = 0;
  u16        * u16tmpimg   = NULL;
  u16        * images      = NULL;
  u8         * results     = NULL;
  rgb_pixel  * rgbtmpimg   = NULL;
  u32          nRows       = 0;
  u32          nCols       = 0;

  /* Warp and registration - k2 lucas kanade */
  fltPixel_t * img        = NULL;
  fltPixel_t * tmplt      = NULL;
  fltPixel_t * IWxp       = NULL;
  fltPixel_t * swap       = NULL;
  u32 M, N;

  int i, a, f, errors;
  u32 nModels, padding_in;

  /* GMM storage */
  u16   * gmm_img    = NULL;
  float * mu         = NULL;
  float * sigma      = NULL;
  float * weight     = NULL;
  u8    * foreground = NULL;

  /* Affine warp parameter set p */
  float warp_p[6];

  if(integrated_binary_read (
    infile, &nRows, &nCols, &padding_in, &M, &N, &IWxp, 
    &nModels, &mu, &sigma, &weight, &nTestImages, &images, &results)) {
    printf("Error! Reading input file failed\n");
    return -1;
  }

  if(padding_in != PAD) {
    printf("Error! Padding in input file does not equal debayer padding...\n");
    return -1;
  }

  STATS_INIT ();

  rgbtmpimg = calloc(M * N, sizeof(rgb_pixel));
  img       = calloc(M * N, sizeof(fltPixel_t));
  tmplt     = calloc(M * N, sizeof(fltPixel_t));

  gmm_img    = calloc(M * N, sizeof(u16));
  foreground = calloc(M * N, sizeof(u8));

  if(!(rgbtmpimg && img && tmplt && IWxp)) {
    fprintf(stderr, "ERROR: Allocating image(s) failed.\n");
  }

  outfp = fopen(outfile, "w");

  for(a = 0; a < nTestImages; a++) {
    PRINT_STAT_INT ("starting_test_image", a);
    /* swap buffers */
    swap = tmplt; tmplt = IWxp; IWxp = swap;

    memset(IWxp, 0, sizeof(fltPixel_t) * M * N);

    /* load next image */
    u16tmpimg = images + nRows * nCols * a;

    /* debay / rgb to luma convert image */
    tic ();
    wami_debayer    (u16tmpimg, rgbtmpimg, nRows, nCols);
    rgb_to_grayscale(rgbtmpimg, img,     M, N);

/* Joseph: Dump grayscale images */
{
   static  int cnt = 0;
   const char name[1024];
   const char namf[1024];
   int i,j ;
   FILE *fp;
   FILE *ff;
   FILE *fr;

   sprintf(name, "./outs_images/grey_img%d.raw", cnt);
   sprintf(namf, "./outs_images/float_img%d.raw", cnt);
   sprintf(namf, "./outs_images/rgb_img%d.raw", cnt);
   cnt ++;
   fp = fopen(namf, "wb+");
   ff = fopen(namf, "wb+");
   fr = fopen(namf, "wb+");
   if((fp == NULL) || (ff == NULL) || (fr == NULL))
   {
      printf("Error in openning the file\n");
      exit(9);
   }
   fwrite(&M, sizeof(M), 1, fp);
   fwrite(&N, sizeof(N), 1, fp);
   fwrite(img, sizeof(*img), M * N, fp);

   fprintf(ff, "%d %d\n", M, N);

   for(i = 0; i < M; ++i){
      for(j = 0; j < N; ++j){
         fprintf(ff, "%g ", img[i*N + j]);
      }
      fprintf(ff, "\n");
   }


   fwrite(&M, sizeof(M), 1, fr);
   fwrite(&N, sizeof(N), 1, fr);
   fwrite(rgbtmpimg, sizeof(*rgbtmpimg), M * N, fr);

   fclose(fp);
   fclose(ff);
   fclose(fr);
}
/* Joseph: Dump grayscale_images */

    PRINT_STAT_DOUBLE ("debayer_and_grayscale", toc ());

print_submatrix_file(outfp, get_name("input_", a), img, M, N, M, N);

    /* register tmplt to previous image */
    for(i = 0; i < 6; i++) warp_p[i] = 0;
    kernel2_lucas_kanade(img, tmplt, warp_p, ITERATIONS, M, N, IWxp);

#ifdef WRITE_OUTPUT_TO_DISK
    print_submatrix_file(outfp, get_name("warp_matrix_", a), warp_p, 2, 3, 2, 3);
#endif

    tic ();
    for(i = 0; i < M * N; i++) { 
      gmm_img[i] = IWxp[i]; 
    }

    f = wami_gmm(M, N, WAMI_GMM_NUM_MODELS, mu, sigma, weight, foreground, gmm_img);
    PRINT_STAT_DOUBLE ("GMM", toc ());
    PRINT_STAT_INT    ("num_foreground", f);

#ifdef WRITE_OUTPUT_TO_DISK
    print_submatrix_u8_file(outfp, get_name("foreground_", a), foreground, M, N, M, N);
#endif

#ifdef ENABLE_CORRECTNESS_CHECKING
    errors = 0;
    for(i = 0; i < M * N; i++) { 
      errors += foreground[i] != results[a * M * N + i];
    }

    PRINT_STAT_INT("num_errors", errors);
    if(((double)errors) > (ERROR_PERCENTAGE * ((double)M * N))) {
      PRINT_STAT_INT("*** SIGNIFICANT ERRORS HAVE OCCURRED IN IMAGE ***", a);
    }
#endif
  }

  STATS_END ();

  free (mu);
  free (sigma);
  free (weight);
  free (foreground);
  free (rgbtmpimg);
  free (images);
  free (results);
  free (img);
  free (tmplt);
  free (IWxp);

  return 0;
}


int
build_test(const char * outfile, int start, int stop, int argc, char * argv[])
{
  /* Input images */
  u16        * u16tmpimg  = NULL;
  rgb_pixel  * rgbtmpimg  = NULL;
  int          nRows      = 0;
  int          nCols      = 0;

  /* Warp and registration - k2 lucas kanade */
  fltPixel_t * img        = NULL;
  fltPixel_t * tmplt      = NULL;
  fltPixel_t * IWxp       = NULL;
  fltPixel_t * swap;
  int M, N;

  int i, a, f;

  /* GMM storage */
  u16   * gmm_img    = NULL;
  float * mu         = NULL;
  float * sigma      = NULL;
  float * weight     = NULL;
  u8    * foreground = NULL;

  /* Affine warp parameter set p */
  float warp_p[6];

  STATS_INIT ();

  /* load first image */
  if(read_octave_u16_2d(argv[1], &nRows, &nCols, &u16tmpimg)) {
    fprintf(stderr, "ERROR: Reading %s failed.\n", argv[1]);
  }

  /* size and allocate remaining workspace */
  M = nRows - 2*PAD;
  N = nCols - 2*PAD;

  rgbtmpimg = calloc(M * N, sizeof(rgb_pixel));
  img       = calloc(M * N, sizeof(fltPixel_t));
  tmplt     = calloc(M * N, sizeof(fltPixel_t));
  IWxp      = calloc(M * N, sizeof(fltPixel_t));

  gmm_img    = calloc(M * N, sizeof(u16));
  mu         = calloc(M * N * WAMI_GMM_NUM_MODELS, sizeof(float));
  sigma      = calloc(M * N * WAMI_GMM_NUM_MODELS, sizeof(float));
  weight     = calloc(M * N * WAMI_GMM_NUM_MODELS, sizeof(float));
  foreground = calloc(M * N, sizeof(u8));

  if(!(rgbtmpimg && img && tmplt && IWxp)) {
    fprintf(stderr, "ERROR: Allocating image(s) failed.\n");
  }

  /* debay / rgb to luma convert first image */
  wami_debayer    (u16tmpimg, rgbtmpimg, nRows, nCols);
  rgb_to_grayscale(rgbtmpimg, IWxp,      M, N);

  for(a = 2; a < argc; a++) {
    PRINT_STAT_INT ("starting_image", a);
    /* swap buffers */
    swap = tmplt; tmplt = IWxp; IWxp = swap;

    memset(IWxp, 0, sizeof(fltPixel_t) * M * N);

    /* load next image */
    tic ();
    if(read_octave_u16_2d(argv[a], &nRows, &nCols, &u16tmpimg)) {
      fprintf(stderr, "ERROR: Reading %s failed.\n", argv[a]);
    }
    PRINT_STAT_DOUBLE ("read_next_image", toc ());

    /* debay / rgb to luma convert first image */
    tic ();
    wami_debayer    (u16tmpimg, rgbtmpimg, nRows, nCols);
    rgb_to_grayscale(rgbtmpimg, img,     M, N);
    PRINT_STAT_DOUBLE ("debayer_and_grayscale", toc ());

    /* register tmplt to previous image */
    for(i = 0; i < 6; i++) warp_p[i] = 0;
    kernel2_lucas_kanade(img, tmplt, warp_p, ITERATIONS, M, N, IWxp);

    tic ();
    for(i = 0; i < M * N; i++) { 
      gmm_img[i] = IWxp[i]; 
    }

    f = wami_gmm(M, N, WAMI_GMM_NUM_MODELS, mu, sigma, weight, foreground, gmm_img);

    PRINT_STAT_DOUBLE ("GMM", toc ());
    PRINT_STAT_INT    ("num_foreground", f);

    if(a == start) {
      PRINT_STAT_INT("Starting test file (image count)", stop - start);
      if(integrated_binary_start (
	  outfile, nRows, nCols, PAD, IWxp, 
	  WAMI_GMM_NUM_MODELS, mu, sigma, weight, 
	  stop - start)) {

	printf("Error creating test file\n");
      }
    } else if (a > start && a <= stop) {
      PRINT_STAT_INT("Writing image", a - start);
      if(integrated_binary_add_image (
	  outfile, nRows, nCols, PAD, u16tmpimg, foreground)) {

	printf("Error writing image to test file\n");
      }
    }
  }

  STATS_END ();

  free (mu);
  free (sigma);
  free (weight);
  free (foreground);
  free (u16tmpimg);
  free (rgbtmpimg);
  free (img);
  free (tmplt);
  free (IWxp);

  return 0;
}
