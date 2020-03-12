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

typedef float fltPixel_t ;

/* compute the 1D x and y gradients of the input matrix / image */
int
gradientXY(fltPixel_t *Iin, int nCols, int nRows, fltPixel_t *Xgradout, fltPixel_t *Ygradout)
{
  int x, y;

  if(!Iin || !Xgradout || !Ygradout)
    return -1;

  if(nCols < 0 || nRows < 0)
    return -2;

  if(!nCols || !nRows)
    return 0;

  /* NOTE this could be made harder to parallelize if the stencil operation iterated over inputs
   * and scattered to output rather than iterating output and accumulating from inputs; however,
   * this could also be made better by collapsing / combining some of these loops */

  /* compute gradient for inside matrix using central difference */
  for (y = 1; y < nRows - 1; y++) {
    for (x = 1; x < nCols - 1; x++) {
      Xgradout[y * nCols + x] =  (  Iin[y     * nCols + (x+1)] - Iin[y     * nCols + (x-1)]) / 2.0;
      Ygradout[y * nCols + x] =  (- Iin[(y-1) * nCols + x    ] + Iin[(y+1) * nCols + x    ]) / 2.0;
    }
  }

  /* handle edge cases */
  /* compute gradient for outer matrix forward / backward difference */
  for (x = 1; x < nCols - 1; x++) {
    Xgradout[x]                   = (Iin[                  (x+1)] - Iin[                  (x-1)]) / 2.0;
    Xgradout[(nRows-1)*nCols + x] = (Iin[(nRows-1)*nCols + (x+1)] - Iin[(nRows-1)*nCols + (x-1)]) / 2.0;
    Ygradout[x]                   = - Iin[x                     ] + Iin[nCols + x];
    Ygradout[(nRows-1)*nCols + x] = Iin[(nRows-1) * nCols + x   ] - Iin[(nRows-2) * nCols + x];
  }
  for (y = 1; y < nRows - 1; y++) {
    Xgradout[y * nCols            ] = - Iin[y * nCols]                + Iin[y * nCols + 1];
    Xgradout[y * nCols + (nCols-1)] = Iin[y * nCols + (nCols-1)]      - Iin[y * nCols + (nCols-2)];
    Ygradout[y * nCols            ] = (Iin[(y+1) * nCols            ] - Iin[(y-1) * nCols            ]) / 2.0;
    Ygradout[y * nCols + (nCols-1)] = (Iin[(y+1) * nCols + (nCols-1)] - Iin[(y-1) * nCols + (nCols-1)]) / 2.0;
  }

  /* compute corners */
  Ygradout [ 0]                             = - Iin [ 0                     ]             + Iin [ nCols + 0];
  Ygradout [ (nRows-1)*nCols + 0]           =   Iin [ (nRows-1) * nCols + 0   ]           - Iin [ (nRows-2) * nCols + 0];
  Ygradout [ (nCols-1)]                     = - Iin [ (nCols-1)                     ]     + Iin [ nCols + (nCols-1)];
  Ygradout [ (nRows-1)*nCols + (nCols-1)]   =   Iin [ (nRows-1) * nCols + (nCols-1)   ]   - Iin [ (nRows-2) * nCols + (nCols-1)];
  Xgradout [ 0 * nCols            ]         = - Iin [ 0 * nCols]                          + Iin [ 0 * nCols + 1];
  Xgradout [ 0 * nCols + (nCols-1)]         =   Iin [ 0 * nCols + (nCols-1)]              - Iin [ 0 * nCols + (nCols-2)];
  Xgradout [ (nRows-1) * nCols            ] = - Iin [ (nRows-1) * nCols]                  + Iin [ (nRows-1) * nCols + 1];
  Xgradout [ (nRows-1) * nCols + (nCols-1)] =   Iin [ (nRows-1) * nCols + (nCols-1)]      - Iin [ (nRows-1) * nCols + (nCols-2)];

  return 0;
}

#ifdef GRADIENT_TEST
int
main(int argc, char ** argv)
{
  float input[] = {0, 1, 0, 1, 0, 1, 0, 1, 0};
  float input_copy[] = {1, 1, 1, 0, 1, 1, 0, 0, 1};
  float inverse[9];
  float result[9];

  printf("Gradient of:\n");
  print_matrix(input, 3, 3);
  gradientXY(input, 3, 3, inverse, result);
  printf("is:\n");
  print_matrix(inverse, 3, 3);
  print_matrix(result, 3, 3);
}
#endif
