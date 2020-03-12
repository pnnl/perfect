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

#include "fileio.h"

/* FILE FORMAT 
 *
 * u32                                  srcRows
 * u32                                  srcCols
 * u32                                  padding
 * u32                                  nRows
 * u32                                  nCols
 * fltPixel_t [nRows x nCols]           start_image (debayed, luma, registered, GMMed)
 * u32                                  nModels
 * float      [nRows x nCols x nModels] mu
 * float      [nRows x nCols x nModels] sigma
 * float      [nRows x nCols x nModels] weight
 * u32                                  nTestImages
 * [nTestImages]
 *   u16      [srcRows x srcCols]       Test Images (bayer filtered)
 *   u8       [srcRows x srcCols]       'correct' result
 */

int
integrated_binary_start (
  const char * filename, 
  u32          srcRows, 
  u32          srcCols, 
  u32          padding, 
  fltPixel_t * start_image,
  u32          nModels,
  float      * mu,
  float      * sigma,
  float      * weight,
  u32          nTestImages)
{
  u32 nRows, nCols;
  
  FILE * fp = fopen(filename, "w");

  if(!fp) {
    return -1;
  }

  fwrite(&srcRows, sizeof(u32), 1, fp);
  fwrite(&srcCols, sizeof(u32), 1, fp);
  fwrite(&padding, sizeof(u32), 1, fp);

  nRows = srcRows - 2*padding;
  nCols = srcCols - 2*padding;

  fwrite(&nRows,       sizeof(u32), 1, fp);
  fwrite(&nCols,       sizeof(u32), 1, fp);
  fwrite(start_image,  sizeof(fltPixel_t), nRows * nCols, fp);
  fwrite(&nModels,     sizeof(u32), 1, fp);
  fwrite(mu,           sizeof(float), nRows * nCols * nModels, fp);
  fwrite(sigma,        sizeof(float), nRows * nCols * nModels, fp);
  fwrite(weight,       sizeof(float), nRows * nCols * nModels, fp);
  fwrite(&nTestImages, sizeof(u32), 1, fp);

  fclose(fp);
  return 0;
}

int
integrated_binary_add_image (
  const char * filename, 
  u32          srcRows, 
  u32          srcCols, 
  u32          padding,
  u16        * image,
  u8         * result)
{
  FILE * fp = fopen(filename, "a");

  if(!fp) {
    return -1;
  }

  fwrite(image,  sizeof(u16), srcRows * srcCols, fp);
  fwrite(result, sizeof(u8),  (srcRows - 2*padding)  * (srcCols - 2*padding), fp);

  fclose(fp);
  return 0;
}

int
integrated_binary_read (
  const char * filename, 
  u32        *  srcRows, 
  u32        *  srcCols, 
  u32        *  padding, 
  u32        *  nRows,
  u32        *  nCols,
  fltPixel_t ** start_image,
  u32        *  nModels,
  float      ** mu,
  float      ** sigma,
  float      ** weight,
  u32        *  nTestImages,
  u16        ** images,
  u8         ** results)
{
  int bytes_read = 0;
  int i = 0;
  FILE * fp = fopen(filename, "r");

  if(!fp) {
    return -1;
  }

  bytes_read += fread(srcRows, sizeof(u32), 1, fp);
  bytes_read += fread(srcCols, sizeof(u32), 1, fp);
  bytes_read += fread(padding, sizeof(u32), 1, fp);
  bytes_read += fread(nRows,   sizeof(u32), 1, fp);
  bytes_read += fread(nCols,   sizeof(u32), 1, fp);

  if(!bytes_read) {
    return -2;
  }
  bytes_read = 0;

  *start_image = malloc(sizeof(fltPixel_t) * *nRows * *nCols);

  if(!start_image) {
    return -1;
  }

  bytes_read += fread(*start_image, sizeof(fltPixel_t), *nRows * *nCols, fp);
  bytes_read += fread(nModels,     sizeof(u32), 1, fp);

  if(!bytes_read) {
    return -2;
  }
  bytes_read = 0;

  *mu     = malloc(sizeof(float) * *nRows * *nCols * *nModels);
  *sigma  = malloc(sizeof(float) * *nRows * *nCols * *nModels);
  *weight = malloc(sizeof(float) * *nRows * *nCols * *nModels);

  if(!mu || ! sigma || !weight) {
    return -1;
  }

  bytes_read += fread(*mu,         sizeof(float), *nRows * *nCols * *nModels, fp);
  bytes_read += fread(*sigma,      sizeof(float), *nRows * *nCols * *nModels, fp);
  bytes_read += fread(*weight,     sizeof(float), *nRows * *nCols * *nModels, fp);
  bytes_read += fread(nTestImages, sizeof(u32), 1, fp);

  if(!bytes_read) {
    return -2;
  }
  bytes_read = 0;

  *images  = malloc(sizeof(u16) * *srcRows * *srcCols * *nTestImages);
  *results = malloc(sizeof(u8)  * *nRows   * *nCols   * *nTestImages);

  if(!images) {
    return -1;
  }

  for(i = 0; i < *nTestImages; i++) {
    bytes_read += fread((*images)  + (i * *srcRows * *srcCols),	
			  sizeof(u16), *srcRows * *srcCols, fp);

    bytes_read += fread((*results) + (i * *nRows * *nCols),  
			  sizeof(u8),  *nRows   * *nCols, fp);
  }

  if(!bytes_read) {
    return -2;
  }

  fclose(fp);

  return 0;
}
