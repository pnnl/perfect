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

#include "sort.h"

int
partition (float * array, int * index, int low, int high)
{
  int left, right, mid;
  int pivot;
  float cur;
  int idx;

  mid = (low + high) / 2;
  left = low;
  right = high;

  /* choose pivot as median of 3: low, high, and mid */
  if ((array[low] - array[mid]) * (array[high] - array[low]) >= 0)
    pivot = low;
  else if ((array[mid] - array[low]) * (array[high] - array[mid]) >= 0)
    pivot = mid;
  else
    pivot = high; 

  /* store value,index at the pivot */
  cur = array[pivot];
  idx = index[pivot];

  /* swap pivot with the first entry in the list */
  array[pivot] = array[low];
  array[low] = cur;
  
  index[pivot] = array[pivot];
  index[low] = idx;

  /* the quicksort itself */
  while (left < right)
  {
    while (array[left] <= cur && left < high)
      left++;
    while (array[right] > cur)
      right--;
    if (left < right)
    {
      float tmp_val;
      int tmp_idx;

      tmp_val = array[right];
      array[right] = array[left];
      array[left] = tmp_val;

      tmp_idx =  index[right];
      index[right] = index[left];
      index[left] = tmp_idx;
    }
  }

  /* pivot was in low, but now moves into position at right */
  array[low] = array[right];
  array[right] = cur;

  index[low] = index[right];
  index[right] = idx;

  return right;
}


/* This defines the length at which we switch to insertion sort */
#define MAX_THRESH 10

int
quicksort_inner (float * array, int * index, int low, int high)
{
  int pivot;
  int length = high - low + 1;

  if (high > low)
  {
    if (length > MAX_THRESH) {
      pivot = partition (array, index, low, high);
      quicksort_inner (array, index, low, pivot-1);
      quicksort_inner (array, index, pivot+1, high);
    }
  }

  return 0;
}

int quicksort (float * array, int * index, int len)
{
  quicksort_inner (array, index, 0, len-1);
  insertion_sort (array, index, len);

  return 0;
}
