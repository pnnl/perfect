/**************************/
/***    UNCLASSIFIED    ***/
/**************************/

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


/***

ALL SOURCE CODE PRESENT IN THIS FILE IS UNCLASSIFIED AND IS
BEING PROVIDED IN SUPPORT OF THE DARPA PERFECT PROGRAM.

THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY, EXPRESSED, IMPLIED, 
OR OTHERWISE INFERRED. USE AND SUITABILITY FOR ANY PARTICULAR
APPLICATION IS SOLELY THE RESPONSIBILITY OF THE IMPLEMENTER. 
NO CLAIM OF SUITABILITY FOR ANY APPLICATION IS MADE.
USE OF THIS CODE FOR ANY APPLICATION RELEASES THE AUTHOR
AND THE US GOVT OF ANY AND ALL LIABILITY.

THE POINT OF CONTACT FOR QUESTIONS REGARDING THIS SOFTWARE IS:

US ARMY RDECOM CERDEC NVESD, RDER-NVS-SI (JOHN HODAPP), 
10221 BURBECK RD, FORT BELVOIR, VA 22060-5806

THIS HEADER SHALL REMAIN PART OF ALL SOURCE CODE FILES.

***/


#include "ensembles.h"


int ensemble_1(algPixel_t *streamA, algPixel_t *out, int nRows, int nCols, int nInpBpp, int nOutBpp)
{	
	// NF
	// H
	// HE 
	// DWT

	int err = 0;
	int nHistBins = 1 << nInpBpp;
	int *h = (int *)calloc(nHistBins, sizeof(int));
	algPixel_t *wrkBuf1 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *wrkBuf2 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	if (!(h && wrkBuf1 && wrkBuf2))
	{
		free(h);
		free(wrkBuf1);
		free(wrkBuf2);
		fprintf(stderr, "File %s, Line %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	memcpy(wrkBuf1, streamA, nRows * nCols * sizeof(algPixel_t));

	err = nf(wrkBuf1, wrkBuf2, nRows, nCols);
	err = hist(wrkBuf2, h, nRows, nCols, nInpBpp);
	err = histEq(wrkBuf2, out, h, nRows, nCols, nInpBpp, nOutBpp);

	err = dwt53(out, nRows, nCols);

	// FOR TESTING, INVERT BACK TO GET DECENT IMAGE FOR COMPARISON...
	// dwt53_inverse(wrkBuf2, nRows, nCols);

	memcpy(out, wrkBuf2, nRows * nCols * sizeof(algPixel_t));
	
	free(wrkBuf2);
	free(wrkBuf1);
	free(h);

	return err;
}


int ensemble_2(algPixel_t *streamA, algPixel_t *streamB, algPixel_t *out, int nRows, int nCols, int nInpBpp, int nOutBpp)
{
	// NF
	// H
	// HE
	// DWT_Fusion

	int err = 0;
	int *histA = (int *)calloc(1 << nInpBpp, sizeof(int));
	int *histB = (int *)calloc(1 << nInpBpp, sizeof(int));

	algPixel_t *wrkBufA = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *wrkBufB = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	if (!(histA && histB && wrkBufA && wrkBufB))
	{
		free(histA);
		free(histB);
		free(wrkBufA);
		free(wrkBufB);
		fprintf(stderr, "File %s, Line %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	err = nf(streamA, wrkBufA, nRows, nCols);
	err = hist(wrkBufA, histA, nRows, nCols, nInpBpp);
	err = histEq(wrkBufA, wrkBufA, histA, nRows, nCols, nInpBpp, nOutBpp);

	err = nf(streamB, wrkBufB, nRows, nCols);
	err = hist(wrkBufB, histB, nRows, nCols, nInpBpp);
	err = histEq(wrkBufB, wrkBufB, histB, nRows, nCols, nInpBpp, nOutBpp);

	err = dwtFusion(wrkBufA, wrkBufB, out, nRows, nCols);

	free(histA);
	free(histB);
	free(wrkBufA);
	free(wrkBufB);

	return err;
}

int ensemble_3(algPixel_t *streamA, algPixel_t *streamB, algPixel_t *out, int nRows, int nCols, int nInpBpp, int nOutBpp)
{
	// NF
	// H
	// HE
	// LAP_Fusion

	int err = 0;
	int *histA = (int *)calloc(1 << nInpBpp, sizeof(int));
	int *histB = (int *)calloc(1 << nInpBpp, sizeof(int));

	algPixel_t *wrkBufA = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *wrkBufB = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	if (!(histA && histB && wrkBufA && wrkBufB))
	{
		free(histA);
		free(histB);
		free(wrkBufA);
		free(wrkBufB);
		fprintf(stderr, "File %s, Line %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	err = nf(streamA, wrkBufA, nRows, nCols);
	err = hist(wrkBufA, histA, nRows, nCols, nInpBpp);
	err = histEq(wrkBufA, wrkBufA, histA, nRows, nCols, nInpBpp, nOutBpp);

	err = nf(streamB, wrkBufB, nRows, nCols);
	err = hist(wrkBufB, histB, nRows, nCols, nInpBpp);
	err = histEq(wrkBufB, wrkBufB, histB, nRows, nCols, nInpBpp, nOutBpp);

#ifdef _USE_DOUBLEDENSITY_LAP_FUSION
	err = lapFusionDD(wrkBufA, wrkBufB, out, nRows, nCols);
#else
	err = lapFusion(wrkBufA, wrkBufB, out, nRows, nCols);
#endif

	free(histA);
	free(histB);
	free(wrkBufA);
	free(wrkBufB);

	return err;
}


