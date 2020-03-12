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

#include "octave_u16.h"

int
read_octave_u16_2d(char * filename, int * nRows, int * nCols, u16 ** output)
{
  FILE * fp = NULL;
  char cur;
  int M, N;
  int r = 0, c = 0;
  int val;
  u16 * out = NULL;
  int read;

  fp = fopen(filename , "r");

  if(!fp) {
    return -1;
  }

  cur = fgetc(fp);
  while(cur == '#' || cur == '%') {
    while(cur != '\n') {
      cur = fgetc(fp);
    }
    cur = fgetc(fp);
  }
  ungetc(cur, fp);

  read = fscanf(fp, "%d %d", &M, &N);

  if(!read) {
    return -1;
  }

  /* image already allocated */
  if(*output) {
    if(M != *nRows || N != *nCols) {
      return -2;
    }
  } else {
    *nRows = M;
    *nCols = N;
  }

  if(!(M * N)) {
    fclose(fp);
    return -3;
  }

  if(!*output) {
    out = malloc(sizeof(u16) * M * N);
  } else {
    out = *output;
  }

  for(c = 0; c < N && !feof(fp); c++) {
    for(r = 0; r < M && !feof(fp); r++) {
      read = fscanf(fp, "%d", &val);
      out[r * N + c] = val;
    }
  }

  if(!read || c < N || r < M) {
    return -4;
  }

  *output = out;
  return 0;
}

int
write_octave_u8_2d(u8 * data, unsigned int nRows, unsigned int nCols, char * filename, char * name)
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
      fprintf(fp, " %d", (int)data[y * nCols + x]);
    }
    fprintf(fp, "\n");
  }

  fclose(fp);

  return 0;
}
