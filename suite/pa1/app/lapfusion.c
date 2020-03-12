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


#include "lapfusion.h"


int lapFusion(algPixel_t *streamA, algPixel_t *streamB, algPixel_t *out, int nRows, int nCols)
{
	int nSubSampRows, nSubSampCols, nSubSampRows2, nSubSampCols2;
	int err = 0;

	fltPixel_t mask[] = {
				 		 1, 1, 1, 1, 1,
						 1, 4, 4, 4, 1,
						 1, 4, 6, 4, 1,
						 1, 4, 4, 4, 1,
						 1, 1, 1, 1, 1
						};


	algPixel_t *imgALowFreq1 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgALowFreq2 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));		
	algPixel_t *imgALowFreq3 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgAHiFreq1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgAHiFreq2  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgAHiFreq3  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *imgBLowFreq1 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgBLowFreq2 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));		
	algPixel_t *imgBLowFreq3 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgBHiFreq1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgBHiFreq2  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgBHiFreq3  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *imgCLowFreq3 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgCLowFreq2 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgCHiFreq3  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgCHiFreq2  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *imgCHiFreq1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *tmp1 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *tmp2 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	if (!(imgALowFreq1 && imgALowFreq2 && imgALowFreq3 &&
		  imgAHiFreq1 && imgAHiFreq2 && imgAHiFreq3 &&
		  imgBLowFreq1 && imgBLowFreq2 && imgBLowFreq3 &&
		  imgBHiFreq1 && imgBHiFreq2 && imgBHiFreq3 &&
		  imgCLowFreq2 && imgCLowFreq3 &&
		  imgCHiFreq1 && imgCHiFreq2 && imgCHiFreq3 &&
		  tmp1 && tmp2))
	{
		// Ok to free() a NULL pointer, so just run through these....
		free(imgALowFreq1); free(imgALowFreq2); free(imgALowFreq3);
		free(imgAHiFreq1); free(imgAHiFreq2); free(imgAHiFreq3);
		free(imgBLowFreq1); free(imgBLowFreq2); free(imgBLowFreq3);
		free(imgBHiFreq1); free(imgBHiFreq2); free(imgBHiFreq3);
		free(imgCLowFreq2); free(imgCLowFreq3);
		free(imgCHiFreq1); free(imgCHiFreq2); free(imgCHiFreq3);
		free(tmp1); 
		free(tmp2);
		fprintf(stderr, "File: %s, Line: %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}


	// Stream A
	err = conv2d(streamA, imgALowFreq1, nRows, nCols, mask, 54.0, 5, 5);
	err = imageSub(streamA, imgALowFreq1, imgAHiFreq1, nRows, nCols);
	err = subsamp(imgALowFreq1, tmp1, nRows, nCols, 2, 2, &nSubSampRows, &nSubSampCols);
	err = conv2d(tmp1, imgALowFreq2, nSubSampRows, nSubSampCols, mask, 54.0, 5, 5);
	err = imageSub(tmp1, imgALowFreq2, imgAHiFreq2, nSubSampRows, nSubSampCols);
	err = subsamp(imgALowFreq2, tmp1, nSubSampRows, nSubSampCols, 2, 2, &nSubSampRows2, &nSubSampCols2);
	err = conv2d(tmp1, imgALowFreq3, nSubSampRows2, nSubSampCols2, mask, 54.0, 5, 5);
	err = imageSub(tmp1, imgALowFreq3, imgAHiFreq3, nSubSampRows2, nSubSampCols2);
	
	// Stream B
	err = conv2d(streamB, imgBLowFreq1, nRows, nCols, mask, 54.0, 5, 5);
	err = imageSub(streamB, imgBLowFreq1, imgBHiFreq1, nRows, nCols);
	err = subsamp(imgBLowFreq1, tmp1, nRows, nCols, 2, 2, &nSubSampRows, &nSubSampCols);
	err = conv2d(tmp1, imgBLowFreq2, nSubSampRows, nSubSampCols, mask, 54.0, 5, 5);
	err = imageSub(tmp1, imgBLowFreq2, imgBHiFreq2, nSubSampRows, nSubSampCols);
	err = subsamp(imgBLowFreq2, tmp1, nSubSampRows, nSubSampCols, 2, 2, &nSubSampRows2, &nSubSampCols2);
	err = conv2d(tmp1, imgBLowFreq3, nSubSampRows2, nSubSampCols2, mask, 54.0, 5, 5);
	err = imageSub(tmp1, imgBLowFreq3, imgBHiFreq3, nSubSampRows2, nSubSampCols2);

	// Fuse...
	err = imageAvg(imgALowFreq3, imgBLowFreq3, imgCLowFreq3, nSubSampRows2, nSubSampCols2);
	err = imageAbsMax(imgAHiFreq3, imgBHiFreq3, imgCHiFreq3, nSubSampRows2, nSubSampCols2);
	err = imageAdd(imgCLowFreq3, imgCHiFreq3, tmp1, nSubSampRows2, nSubSampCols2);
	err = interp2(tmp1, tmp2, nSubSampRows2, nSubSampCols2, &nSubSampRows, &nSubSampCols);
	err = imageAbsMax(imgAHiFreq2,imgBHiFreq2, imgCHiFreq2, nSubSampRows, nSubSampCols);
	err = imageAdd(imgCHiFreq2, tmp2, tmp1, nSubSampRows, nSubSampCols);
	err = interp2(tmp1, tmp2, nSubSampRows, nSubSampCols, &nRows, &nCols);
	err = imageAbsMax(imgAHiFreq1, imgBHiFreq1, imgCHiFreq1, nRows, nCols);
	err = imageAdd(tmp2, imgCHiFreq1, out, nRows, nCols);

	// Magic fused answer in 'out[]'

	// Free up all the temp buffers we used...

	free(imgALowFreq1);
	free(imgALowFreq2);
	free(imgALowFreq3);
	free(imgAHiFreq1);
	free(imgAHiFreq2);
	free(imgAHiFreq3);

	free(imgBLowFreq1);
	free(imgBLowFreq2);
	free(imgBLowFreq3);
	free(imgBHiFreq1);
	free(imgBHiFreq2);
	free(imgBHiFreq3);

	free(imgCLowFreq3);
	free(imgCLowFreq2);
	free(imgCHiFreq3);
	free(imgCHiFreq2);
	free(imgCHiFreq1);

	free(tmp1);
	free(tmp2);

	return 0;
}


int lapFusionDD(algPixel_t *streamA, algPixel_t *streamB, algPixel_t *out, int nRows, int nCols)
{
	int nSubSampRows, nSubSampCols, nSubSampRows2, nSubSampCols2;
	int nFilterRowsF = 5; 
	int nFilterColsF = 5;
	int err = 0;
	int i = 0;

	fltPixel_t F[] =    {
				 	     1,   4,   6,   4,  1,
	                     4,  16,  24,  16,  4,
                         6,  24,  36,  24,  6,
	                     4,  16,  24,  16,  4,
                         1,   4,   6,   4,  1
						};

	for (i = 0; i < nFilterRowsF * nFilterColsF; i++)
	{
		F[i] /= 256.0;
	}

	int nFilterRowsFD = 9; 
	int nFilterColsFD = 9;
	
	fltPixel_t FD[] =  {
						 1,   3,   4,   5,   6,   5,  4,    3,  1,
						 3,   9,  12,  15,  18,  15,  12,   9,  3,
						 4,  12,  16,  20,  24,  20,  16,  12,  4,
						 5,  15,  20,  25,  30,  25,  20,  15,  5,
						 6,  18,  24,  30,  36,  30,  24,  18,  6,
						 5,  15,  20,  25,  30,  25,  20,  15,  5,
						 4,  12,  16,  20,  24,  20,  16,  12,  4,
						 3,   9,  12,  15,  18,  15,  12,   9,  3,
						 1,   3,   4,   5,   6,   5,   4,   3,  1
						};

	for (i = 0; i < nFilterRowsFD * nFilterColsFD; i++)
	{
		FD[i] /= (1024.0);
	}


	algPixel_t *LF_A1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *LF_A2D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));		
	algPixel_t *LF_A3D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *LF_A4D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *HF_A1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *HF_A2D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));		
	algPixel_t *HF_A3D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *HF_A4D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *LF_B1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *LF_B2D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));		
	algPixel_t *LF_B3D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *LF_B4D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *HF_B1  = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *HF_B2D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));		
	algPixel_t *HF_B3D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *HF_B4D = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	algPixel_t *tmp1 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *tmp2 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *tmp3 = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));

	if (!(LF_A1 && LF_A2D && LF_A3D && LF_A4D &&
		  HF_A1 && HF_A2D && HF_A3D && HF_A4D &&
		  LF_B1 && LF_B2D && LF_B3D && LF_B4D &&
		  HF_B1 && HF_B2D && HF_B3D && HF_B4D &&
		  tmp1 && tmp2 && tmp3))
	{
		free(LF_A1); free(LF_A2D); free(LF_A3D); free(LF_A4D);
		free(HF_A1); free(HF_A2D); free(HF_A3D); free(HF_A4D);
		free(LF_B1); free(LF_B2D); free(LF_B3D); free(LF_B4D);
		free(HF_B1); free(HF_B2D); free(HF_B3D); free(HF_B4D);
		free(tmp1);
		free(tmp2);
		free(tmp3);
		fprintf(stderr, "File: %s, Line: %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	// Stream A
	// LF_A1 = conv2(imA, f, 'same');
	// HF_A1 = imA - LF_A1;
	err = conv2d(streamA, LF_A1, nRows, nCols, F, 1.0, nFilterRowsF, nFilterColsF);
	err = imageSub(streamA, LF_A1, HF_A1, nRows, nCols);

	// LF_A2D = conv2(LF_A1, fd, 'same');
	// HF_A2D = LF_A1 - LF_A2D;
	err = conv2d(LF_A1, LF_A2D, nRows, nCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageSub(LF_A1, LF_A2D, HF_A2D, nRows, nCols);

	// tmp = LF_A2D(1:2:end, 1:2:end);
	err = subsamp(LF_A2D, tmp1, nRows, nCols, 2, 2, &nSubSampRows, &nSubSampCols);

	// LF_A3D   = conv2(tmp, fd, 'same');
	// HF_A3D   = tmp - LF_A3D;
	err = conv2d(tmp1, LF_A3D, nSubSampRows, nSubSampCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageSub(tmp1, LF_A3D, HF_A3D, nSubSampRows, nSubSampCols);

	// tmp = LF_A3D(1:2:end, 1:2:end);
	err = subsamp(LF_A3D, tmp2, nSubSampRows, nSubSampCols, 2, 2, &nSubSampRows2, &nSubSampCols2);

    //LF_A4D   = conv2(tmp, fd, 'same');
	//HF_A4D   = tmp - LF_A4D;
	err = conv2d(tmp2, LF_A4D, nSubSampRows2, nSubSampCols2, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageSub(tmp2, LF_A4D, HF_A4D, nSubSampRows2, nSubSampCols2);

	
	// Stream B

	// LF_B1 = conv2(imB, f, 'same');
	// HF_B1 = imB - LF_B1;
	err = conv2d(streamB, LF_B1, nRows, nCols, F, 1.0, nFilterRowsF, nFilterColsF);
	err = imageSub(streamB, LF_B1, HF_B1, nRows, nCols);

	// LF_B2D = conv2(LF_B1, fd, 'same');
	// HF_B2D = LF_B1 - LF_B2D;
	err = conv2d(LF_B1, LF_B2D, nRows, nCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageSub(LF_B1, LF_B2D, HF_B2D, nRows, nCols);

	// tmp = LF_B2D(1:2:end, 1:2:end);
	err = subsamp(LF_B2D, tmp1, nRows, nCols, 2, 2, &nSubSampRows, &nSubSampCols);

	// LF_B3D   = conv2(tmp, fd, 'same');
	// HF_B3D   = tmp - LF_B3D;
	err = conv2d(tmp1, LF_B3D, nSubSampRows, nSubSampCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageSub(tmp1, LF_B3D, HF_B3D, nSubSampRows, nSubSampCols);

	// tmp = LF_B3D(1:2:end, 1:2:end);
	err = subsamp(LF_B3D, tmp2, nSubSampRows, nSubSampCols, 2, 2, &nSubSampRows2, &nSubSampCols2);
	
    //LF_B4D   = conv2(tmp, fd, 'same');
	//HF_B4D   = tmp - LF_B4D;
	err = conv2d(tmp2, LF_B4D, nSubSampRows2, nSubSampCols2, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageSub(tmp2, LF_B4D, HF_B4D, nSubSampRows2, nSubSampCols2);


	// Fuse...
    // b	= (LF_A4D + LF_B4D) / 2;
    // b	= b + absMax(HF_A4D, HF_B4D);
	err = imageAvg(LF_A4D, LF_B4D, tmp1, nSubSampRows2, nSubSampCols2);
	err = imageAbsMax(HF_A4D, HF_B4D, tmp2, nSubSampRows2, nSubSampCols2);
	err = imageAdd(tmp1, tmp2, tmp3, nSubSampRows2, nSubSampCols2);

	// upsample...
    // different interpolation function. See function for details.
    // bb	= interpolator(b, fd);
    // bb	= bb + absMax(HF_A3D, HF_B3D) + ...
    //             conv2(absMax(HF_A3D, HF_B3D), fd, 'same');

	err = interpolator(tmp3, tmp1, nSubSampRows2, nSubSampCols2, FD, 1, nFilterRowsFD, nFilterColsFD);
	err = imageAbsMax(HF_A3D, HF_B3D, tmp2, nSubSampRows, nSubSampCols);
	err = conv2d(tmp2, tmp3, nSubSampRows, nSubSampCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageAdd(tmp1, tmp3, tmp2, nSubSampRows, nSubSampCols);


    // upsample...
    // bbb  = interpolator(bb, fd);
    // bbb	= bbb + absMax(HF_A2D, HF_B2D) + ...
    //        conv2(absMax(HF_A2D, HF_B2D), fd, 'same');

	err = interpolator(tmp2, tmp1, nSubSampRows, nSubSampCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageAbsMax(HF_A2D, HF_B2D, tmp2, nRows, nCols);
	err = conv2d(tmp2, tmp3, nRows, nCols, FD, 1.0, nFilterRowsFD, nFilterColsFD);
	err = imageAdd(tmp1, tmp3, tmp2, nRows, nCols);


    // upsample...
    // No upsampling here since the "2nd DD level" was not downsampled.
	// fused = bbb + absMax(HF_A1, HF_B1);
	err = imageAbsMax(HF_A1, HF_B1, tmp1, nRows, nCols);
	err = imageAdd(tmp1, tmp2, out, nRows, nCols);

	// Magic fused answer in 'out[]'

	// Free up all the temp buffers we used...

	free(LF_A1);
	free(LF_A2D);
	free(LF_A3D);
	free(LF_A4D);
	free(LF_B1);
	free(LF_B2D);
	free(LF_B3D);
	free(LF_B4D);

	free(HF_A1);
	free(HF_A2D);
	free(HF_A3D);
	free(HF_A4D);
	free(HF_B1);
	free(HF_B2D);
	free(HF_B3D);
	free(HF_B4D);

	free(tmp1);
	free(tmp2);
	free(tmp3);

	return err;
}
