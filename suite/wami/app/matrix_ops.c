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
#include <math.h>
#include "matrix_ops.h"

#define CLOSE_TO_ZERO(X) (((X) < 1e-36) && ((X) > -1e-36))

int
reshape(fltPixel_t * in, int nRows, int nCols, int newRows, int newCols, fltPixel_t * out)
{
  int cur_index;

  if(!in || !out)
    return -1;

  if(nCols < 0 || nRows < 0 || newRows < 0 || newCols < 0)
    return -2;

  if(newRows * newCols > nRows * nCols)
    return -3;

  if(!nCols || !nRows || !nRows || !newCols)
    return 0;

  for (cur_index = 0; cur_index < newRows * newCols; cur_index++) {
    out[(cur_index % newRows) * newCols + (cur_index / newRows)] =  in[(cur_index % nRows) * nCols + (cur_index / nRows)];
  }

  return 0;
}

int
subtract(fltPixel_t * inA, fltPixel_t * inB, int nRows, int nCols, fltPixel_t * out)
{
  int x, y;

  if(!inA|| !inB|| !out)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      out[y * nCols + x] =  inA[y * nCols + x] - inB[y * nCols + x];
    }
  }

  return 0;
}

int
add(fltPixel_t * inA, fltPixel_t * inB, int nRows, int nCols, fltPixel_t * out)
{
  int x, y;

  if(!inA|| !inB|| !out)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      out[y * nCols + x] =  inA[y * nCols + x] + inB[y * nCols + x];
    }
  }

  return 0;
}

/* A [nRows x nCommon] * B [nCommon x nCols] -> out[nRows x nCols] */
int
mult(fltPixel_t * inA, fltPixel_t * inB, int nRows, int nCols, int nCommon, fltPixel_t * out)
{
  int x, y, i;

  if(!inA|| !inB|| !out)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      out[y * nCols + x] = 0;
    }
  }

  /* most trivial implementation */
  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      for (i = 0; i < nCommon; i++) {
        out[y * nCols + x] += inA[y * nCommon + i] * inB[i * nCols + x];
      }
    }
  }

  return 0;
}

int
swap_row(fltPixel_t * data, int nRows, int nCols, int r1, int r2)
{
  int x;

  if(!data)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(r1 > nRows || r2 > nRows)
    return -3;

  if(!nCols || !nRows)
    return 0;

  for (x = 0; x < nCols; x++) {
    float tmp = data[r1 * nCols + x];
    data[r1 * nCols + x] = data[r2 * nCols + x];
    data[r2 * nCols + x] = tmp;
  }

  return 0;
}

int
scale_row(fltPixel_t * data, int nRows, int nCols, int r, float scale)
{
  int x;

  if(!data)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(r > nRows)
    return -3;

  if(!nCols || !nRows)
    return 0;

  for (x = 0; x < nCols; x++) {
    data[r * nCols + x] *= scale;
  }

  return 0;
}

int
scale_and_add_row(fltPixel_t * data, int nRows, int nCols, int r1, int r2, float scale)
{
  int x;

  if(!data)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(r1 > nRows || r2 > nRows)
    return -3;

  if(!nCols || !nRows)
    return 0;

  for (x = 0; x < nCols; x++) {
    data[r1 * nCols + x] += scale * data[r2 * nCols + x];
  }

  return 0;
}

int
invert_gauss_jordan(fltPixel_t * data, fltPixel_t * workspace, int nRows, int nCols, fltPixel_t * inverse)
{
  int x, y;
  int r, c;

  if(!data || !inverse)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;
 
  /* set inverse to identity */
  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      inverse  [y * nCols + x] = (x == y) ? 1.0 : 0.0;
      workspace[y * nCols + x] = data[y * nCols + x];
    }
  }

  /* for each column */
  for (c = 0; c < nCols; c++) {
    float scale = 1.0;

    /* swap rows if close to zero */
    if(CLOSE_TO_ZERO(workspace[c * nCols + c])) {
      for(r = c + 1; r < nCols; r++) {
	if(!CLOSE_TO_ZERO(workspace[r * nCols + c])) {
	  swap_row(workspace, nRows, nCols, r, c);
	  swap_row(inverse, nRows, nCols, r, c);
	  break;
	}
      }
      if(r >= nCols)
	return -3;
     }

    /* scale operation */
    scale = 1.0f / workspace[c * nCols + c];
    scale_row(workspace, nRows, nCols, c, scale);
    scale_row(inverse, nRows, nCols, c, scale);

    scale = 1.0f / workspace[c * nCols + c];
    scale_row(workspace, nRows, nCols, c, scale);
    scale_row(inverse, nRows, nCols, c, scale);

    /* zero column */
    for(r = 0; r < nRows; r++) {
      if(r != c) {
        scale = - workspace[r * nCols + c];
        scale_and_add_row(workspace, nRows, nCols, r, c, scale);
        scale_and_add_row(inverse, nRows, nCols, r, c, scale);
      }
    }
  }

  return 0;
}

int
print_matrix(fltPixel_t * mat, int nRows, int nCols)
{
  int x, y;

  if(!mat)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  printf("[");
  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      printf(" %8e ", mat[y * nCols + x]);
    }
    printf(y == (nRows-1) ? "]\n" : ";\n");
  }

  return 0;
}

int
print_submatrix(fltPixel_t * mat, int nRows, int nCols, int width, int height)
{
  int x, y;

  if(!mat)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  printf("[\n");
  for (y = 0; y < width; y++) {
    for (x = 0; x < height; x++) {
      printf(" %8e ", mat[y * nCols + x]);
    }
    printf(";\n");
  }
  printf("]\n");

  return 0;
}

int
print_submatrix_u16(u16 * mat, int nRows, int nCols, int width, int height)
{
  int x, y;

  if(!mat)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  printf("[\n");
  for (y = 0; y < width; y++) {
    for (x = 0; x < height; x++) {
      printf(" %d ", (int)mat[y * nCols + x]);
    }
    printf(";\n");
  }
  printf("]\n");

  return 0;
}

int
print_submatrix_u8_file(FILE * fp, char * name, u8 * mat, int nRows, int nCols, int width, int height)
{
  int x, y;

  if(!mat || !fp)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  if (name) {
    fprintf(fp,"%s = [\n", name);
  } else {
    fprintf(fp,"[\n");
  }
  for (y = 0; y < width; y++) {
    for (x = 0; x < height; x++) {
      fprintf(fp," %d ", (int)mat[y * nCols + x]);
    }
    fprintf(fp, ";\n");
  }
  fprintf(fp, "];\n");

  return 0;
}

int
print_submatrix_file(FILE * fp, char * name, fltPixel_t * mat, int nRows, int nCols, int width, int height)
{
  int x, y;

  if(!mat || !fp)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  if (name) {
    fprintf(fp,"%s = [\n", name);
  } else {
    fprintf(fp,"[\n");
  }
  for (y = 0; y < width; y++) {
    for (x = 0; x < height; x++) {
      fprintf(fp," %8e ", (float)mat[y * nCols + x]);
    }
    fprintf(fp, ";\n");
  }
  fprintf(fp, "];\n");

  return 0;
}

int
print_submatrix_rgb(rgb_pixel * mat, int nRows, int nCols, int width, int height, int channel)
{
  int x, y;

  if(!mat)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  printf("[\n");
  switch(channel) {
    case 0:
      for (y = 0; y < width; y++) {
	for (x = 0; x < height; x++) {
	  printf(" %d ", (int)mat[y * nCols + x].r);
	}
	printf(";\n");
      }
    break;

    default:
    case 1:
      for (y = 0; y < width; y++) {
	for (x = 0; x < height; x++) {
	  printf(" %d ", (int)mat[y * nCols + x].g);
	}
	printf(";\n");
      }
    break;

    case 2:
      for (y = 0; y < width; y++) {
	for (x = 0; x < height; x++) {
	  printf(" %d ", (int)mat[y * nCols + x].b);
	}
	printf(";\n");
      }
    break;
  }
  printf("]\n");

  return 0;
}

#ifdef MATRIX_OPS_TEST
int
main(int argc, char ** argv)
{
  float input[] = {1, 1, 1, 0, 1, 1, 0, 0, 1};
  float input_copy[] = {1, 1, 1, 0, 1, 1, 0, 0, 1};
  float inverse[9];
  float result[9];

  printf("Inverse of:\n");
  print_matrix(input, 3, 3);
  invert_gauss_jordan(input_copy, 3, 3, inverse);
  printf("is:\n");
  print_matrix(inverse, 3, 3);
  mult(input, inverse, 3, 3, 3, result);
  printf("multiplied:\n");
  print_matrix(result, 3, 3);
  add(input, inverse, 3, 3, result);
  printf("added:\n");
  print_matrix(result, 3, 3);
  subtract(input, inverse, 3, 3, result);
  printf("subtracted:\n");
  print_matrix(result, 3, 3);

  float mata[] = {1, 1, 1, 0, 1, 1, 0, 0, 1};
  float matb[] = {1, 1, 1, 0, 1, 1};
  float matc[] = {0,0,0,0,0,0};

  printf("mata is:\n");
  print_matrix(mata, 3,3);
  printf("matb is:\n");
  print_matrix(matb, 3,2);
  mult(mata, matb, 3,2,3,matc);
  printf("matc is:\n");
  print_matrix(matc, 3,2);
}
#endif
