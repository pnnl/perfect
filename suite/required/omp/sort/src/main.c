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
#include <time.h>
#include <omp.h>

#include "sort.h"
#include "timing/timer.h"
#include "xmem/xmalloc.h"

#if !defined(TAV_SORT_BATCH_SIZE)
#define TAV_SORT_BATCH_SIZE (1024)
#endif

#if INPUT_SIZE == INPUT_SIZE_SMALL
#define TAV_SORT_LENGTH (9)

#elif INPUT_SIZE == INPUT_SIZE_MEDIUM
#define TAV_SORT_LENGTH (1<<6)

#elif INPUT_SIZE == INPUT_SIZE_LARGE
#define TAV_SORT_LENGTH (1<<10) 

#else
#error "Unhandled value for INPUT_SIZE"
#endif


/* #define USE_N_SQUARED_SORT (1) */
/* #define USE_RADIX_SORT (1) */

int main (int argc, char ** argv)
{
  double t_start, t_end;
  int i, j;
  int * array_index;
  float * array_value;
  int err;

  srand (time (NULL));

  /* Print parameters */
  STATS_INIT ();
  PRINT_STAT_STRING ("kernel", "sort");
  PRINT_STAT_INT ("length", TAV_SORT_LENGTH);
  PRINT_STAT_INT ("batch_size", TAV_SORT_BATCH_SIZE);

  /* Allocate TAV_SORT_LENGTH arrays for values (to sort by) and indices (metadata) */
  tic ();
  array_index = (int *) xmalloc (sizeof(int) * TAV_SORT_BATCH_SIZE * TAV_SORT_LENGTH);
  array_value = (float *) xmalloc (sizeof(float) * TAV_SORT_BATCH_SIZE * TAV_SORT_LENGTH);
  PRINT_STAT_DOUBLE ("time_malloc", toc ());

#if defined(USE_RADIX_SORT)
  int * fake_array;
  fake_array = (int *) xmalloc (sizeof(int) * TAV_SORT_BATCH_SIZE * TAV_SORT_LENGTH);
#endif

  /* Initialize with random values */
  tic ();
  for (j = 0; j < TAV_SORT_BATCH_SIZE; j++)
  {
    for (i = 0; i < TAV_SORT_LENGTH; i++)
    {
      array_index[j*TAV_SORT_LENGTH + i] = i;
      array_value[j*TAV_SORT_LENGTH + i] = ((float) rand () / (float) RAND_MAX) ;
#if defined(USE_RADIX_SORT)
    fake_array[j*TAV_SORT_LENGTH + i] =  (int)(10000 *  array_value[j*TAV_SORT_LENGTH + i]);
#endif
    }
  }
  PRINT_STAT_DOUBLE ("time_generate_random_data", toc ());

  /* Do the actual sort -- several implementations provided */
  tic ();

  for (j = 0; j < TAV_SORT_BATCH_SIZE; j++)
  {
    float * array_val_to_sort = &array_value[j*TAV_SORT_LENGTH];
    int * array_idx_to_sort = &array_index[j*TAV_SORT_LENGTH];
#if defined(USE_N_SQUARED_SORT)
    err = n_squared_sort (array_val_to_sort, array_idx_to_sort, TAV_SORT_LENGTH);
#elif defined(USE_RADIX_SORT)
    err = radix_sort_tuples (&fake_array[j*TAV_SORT_LENGTH], array_idx_to_sort, TAV_SORT_LENGTH, 8);
#elif defined(USE_INSERTION_SORT)
    err = insertion_sort (array_val_to_sort, array_idx_to_sort, TAV_SORT_LENGTH);
#else
    err = quicksort (array_val_to_sort, array_idx_to_sort, TAV_SORT_LENGTH);
#endif
  }
  PRINT_STAT_DOUBLE ("time_sort", toc ());

  /* Validate the sorting -- should return 0 if all arrays sorted correctly */
  tic ();
  err = 0;
  for (j = 0; j < TAV_SORT_BATCH_SIZE; j++)
  {
    float * array_val_sorted = &array_value[j*TAV_SORT_LENGTH];
#if defined(USE_RADIX_SORT)
    array_val_sorted = &fake_array[j*TAV_SORT_LENGTH];
#endif
    err += validate_sorted (array_val_sorted, TAV_SORT_LENGTH);
  }
  PRINT_STAT_DOUBLE ("time_validate", toc ());
  PRINT_STAT_INT ("result_validate (number of errors)", err);

  STATS_END ();

  free (array_value);
  free (array_index);

  return 0;
}
