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
#include "octave.h"

int
write_array_to_octave (int * data, int nrows, int ncols, char * filename, char * name)
{
  FILE *fp = fopen(filename, "w");
  int i, j;
  int n = nrows;
  int m = ncols;

  fprintf(fp, "%% Created by PERFECT Benchmark Suite\n");
  fprintf(fp, "%% name: %s\n", name);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %d\n", n);
  fprintf(fp, "%% columns: %d\n", m);

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      fprintf(fp, " %d", data[i*m+j]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  return 0;
}

int
read_array_from_octave (int * data, int nrows, int ncols, char * filename)
{
  FILE *fp = fopen(filename, "r");
  int i, j;
  int n = nrows;
  int m = ncols;
  char buffer[100];
  char comment;
  int rc= 0;
  char *rc_ptr = NULL;

  if (!fp) {
    fprintf(stderr, "File not found: %s\n", filename);
    exit(-1);
  }

  rc = fscanf(fp, "%c", &comment);
  while (comment == '#' || comment == '%') {
    rc_ptr = fgets(buffer, 100, fp);
    rc = fscanf(fp, "%c", &comment);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      int tmp;
      rc = fscanf(fp, "%d", &tmp);
      data[i*m+j] = (int) tmp;
    }
    rc = fscanf(fp, "\n");
  }

  fclose(fp);

  return 0;
}

int
write_fltarray_to_octave (float * data, unsigned int nRows, unsigned int nCols, char * filename, char * name)
{
  FILE *fp = fopen(filename, "w");
  int x, y;

  fprintf(fp, "%% Created by PERFECT Benchmark Suite\n");
  fprintf(fp, "%% name: %s\n", name);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %d\n", nRows);
  fprintf(fp, "%% columns: %d\n", nCols);

  for (y = 0; y < nRows; y++) {
    for (x = 0; x < nCols; x++) {
      fprintf(fp, " %.7g", data[y * nCols + x]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  return 0;
}

int
read_fltarray_from_octave (float * data, int nrows, int ncols, char * filename)
{
  FILE *fp = fopen(filename, "r");
  int i, j;
  int n = nrows;
  int m = ncols;
  char buffer[100];
  char comment;
  int rc = 0;
  char *rc_ptr = NULL;

  if (!fp) {
    fprintf(stderr, "File not found: %s\n", filename);
    exit(-1);
  }

  rc = fscanf(fp, "%c", &comment);
  while (comment == '#' || comment == '%') {
    rc_ptr = fgets(buffer, 100, fp);
    rc = fscanf(fp, "%c", &comment);
  }

  for (i = 0; i < n; i++) {
    for (j = 0; j < m; j++) {
      float tmp;
      rc = fscanf(fp, "%f", &tmp);
      data[i*m+j] = (float) tmp;
    }
    rc = fscanf(fp, "\n");
  }

  fclose(fp);

  return 0;
}

int
write_rgbarray_to_octave (rgb_pixel * data, unsigned int nrows, unsigned int ncols, int channel, char * filename, char * name)
{
  FILE *fp = fopen(filename, "w");
  int x, y;

  fprintf(fp, "%% Created by PERFECT Benchmark Suite\n");
  fprintf(fp, "%% name: %s\n", name);
  fprintf(fp, "%% type: matrix\n");
  fprintf(fp, "%% rows: %d\n", nrows);
  fprintf(fp, "%% columns: %d\n", ncols);

  switch(channel) {
    case 0:
      for (y = 0; y < nrows; y++) {
	for (x = 0; x < ncols; x++) {
	  fprintf(fp, " %d", ((int)data[y * ncols + x].r));
	}
	fprintf(fp, "\n");
      }
    break;

    default:
    case 1:
      for (y = 0; y < nrows; y++) {
	for (x = 0; x < ncols; x++) {
	  fprintf(fp, " %d", ((int)data[y * ncols + x].g));
	}
	fprintf(fp, "\n");
      }
    break;

    case 2:
      for (y = 0; y < nrows; y++) {
	for (x = 0; x < ncols; x++) {
	  fprintf(fp, " %d", ((int)data[y * ncols + x].b));
	}
	fprintf(fp, "\n");
      }
    break;
  }

  fclose(fp);

  return 0;
}
