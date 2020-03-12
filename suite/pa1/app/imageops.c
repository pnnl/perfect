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


#include "imageops.h"


int imageStats(algPixel_t *image, float *imgMean, float *imgStd, float *imgMax, float *imgMin, int nRows, int nCols)
{
	int nPxls = nRows * nCols;
	int i = 0;
	float diffSq = 0.0;
	float mean = 0.0; 
	float stdv = 0.0;
	float sum = 0.0;
	float min = (float)INT_MAX;
	float max = (float)INT_MIN;

	if (nPxls <= 1)
	{
		return -1;	// Bad input
	}

	if (nPxls <= 1)
	{
		*imgMean = 0;
		*imgStd = 0;
		return -1;
	}

	for (i = 0; i < nPxls; i++)
	{
		sum += (float)image[i];
	}

	mean = sum / nPxls;

	for (i = 0; i < nPxls; i++)
	{
		diffSq += ((image[i] - mean) * (image[i] - mean));
		min = MIN(min, image[i]);
		max = MAX(max, image[i]);
	}
	stdv = sqrt((diffSq / (nPxls - 1)));

	*imgMean = mean;
	*imgStd  = stdv;
	*imgMin  = min;
	*imgMax  = max;
 
	return 0;
}


int rescaleImage(algPixel_t *image, int nRows, int nCols, int nBpp)
{
	int nGrayLevels = (1 << nBpp) - 1;
	int nPxls = nRows * nCols;
	float fMin = 0.0, fMax = 0.0;
	float fMean = 0.0, fStd = 0.0;
	int err = 0;
	int i = 0;

	err = imageStats(image, &fMean, &fStd, &fMax, &fMin, nRows, nCols);

	for (i = 0; i < nPxls; i++)
	{
		image[i] = (algPixel_t) ((((float)image[i] - fMin) / (fMax - fMin)) * nGrayLevels);
	}

	return err;
}


int imageAdd(algPixel_t *imageA, algPixel_t *imageB, algPixel_t *imageC, int nRows, int nCols)
{
	algPixel_t *pA = imageA;
	algPixel_t *pB = imageB;
	algPixel_t *pC = imageC;
	int i = 0;

	for (i = 0; i < nRows * nCols; i++)
	{
		*pC++ = (*pA++) + (*pB++);
	}

	return 0;
}

int imageSub(algPixel_t *imageA, algPixel_t *imageB, algPixel_t *imageC, int nRows, int nCols)
{
	algPixel_t *pA = imageA;
	algPixel_t *pB = imageB;
	algPixel_t *pC = imageC;
	int i = 0;

	for (i = 0; i < nRows * nCols; i++)
	{
		*pC++ = (*pA++) - (*pB++);
	}

	return 0;
}

int imageAvg(algPixel_t *imageA, algPixel_t *imageB, algPixel_t *imageC, int nRows, int nCols)
{
	algPixel_t *pA = imageA;
	algPixel_t *pB = imageB;
	algPixel_t *pC = imageC;
	int i = 0;

	for (i = 0; i < nRows * nCols; i++)
	{
		*pC++ = ((*pA++) + (*pB++)) >> 1;
	}

	return 0;
}


int imageMax(algPixel_t *imageA, algPixel_t *imageB, algPixel_t *imageC, int nRows, int nCols)
{
	int i = 0;

	for (i = 0; i < nRows * nCols; i++)
	{
		imageC[i] = MAX(imageA[i], imageB[i]);
	}

	return 0;
}

int imageAbsMax(algPixel_t *imageA, algPixel_t *imageB, algPixel_t *imageC, int nRows, int nCols)
{
	int i = 0;

	for (i = 0; i < nRows * nCols; i++)
	{
		if (ABS(imageA[i]) > ABS(imageB[i]))
		{
			imageC[i] = imageA[i];
		}
		else
		{
			imageC[i] = imageB[i];
		}
	}

	return 0;
}


int imageMat2Gray(algPixel_t *srcImage, algPixel_t *dstImage, int nBpp, int nRows, int nCols)
{
	int err = 0;
	int i = 0;
	//fltPixel_t fPxl;
	algPixel_t *pSrc = srcImage;
	algPixel_t *pDst = dstImage;
	algPixel_t *pTmp = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
	algPixel_t *pBuf = pTmp;

	if (!pTmp)
	{
		fprintf(stderr, "File: %s, Line: %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	int nGrayLevels = (1 << nBpp) - 1; // was (1 << nBpp)
	int nPxls = nRows * nCols;
	int low = 0, high = 0, delta = 0, Pxl = 0;
	float mean = 0.0, stdv = 0.0;
	float min = 0.0, max = 0.0;

	err = imageStats(srcImage, &mean, &stdv, &max, &min, nRows, nCols);

	low  = (int) (mean - 2 * stdv - 0.5);
	if (low < 0)
	{
		low = 0;
	}

	high = (int) (mean + 2 * stdv + 0.5);
	delta = (high - low);
	if (delta == 0)
	{
		delta = 1;
	}

	for (i = 0; i < nPxls; i++)
	{
		Pxl = *pSrc++;
		if (Pxl < low)
		{
			Pxl = low;
		}
		else if (Pxl > high)
		{
			Pxl = high;
		}

		*pBuf++ = (algPixel_t) (((Pxl - low) / (delta)) * nGrayLevels);
	}

	memcpy(pDst, pTmp, nRows * nCols * sizeof(algPixel_t));
	free(pTmp);

	return err;
}


//int imageMat2Gray(algPixel_t *srcImage, algPixel_t *dstImage, int nBpp, int nRows, int nCols)
//{
//	int err = 0;
//	fltPixel_t fPxl;
//	algPixel_t *pSrc = srcImage;
//	algPixel_t *pDst = dstImage;
//	algPixel_t *pTmp = (algPixel_t *)calloc(nRows * nCols, sizeof(algPixel_t));
//	algPixel_t *pBuf = pTmp;
//
//	if (!pTmp)
//	{
//		fprintf(stderr, "File: %s, Line: %d, Memory Allocation Error\n", __FILE__, __LINE__);
//		return -1;
//	}
//
//	int nGrayLevels = (1 << nBpp) - 1; // was (1 << nBpp)
//	int nPxls = nRows * nCols;
//	float low = 0.0, high = 0.0;
//	float mean = 0.0, stdv = 0.0;
//	float min = 0.0, max = 0.0;
//
//	err = imageStats(srcImage, &mean, &stdv, &max, &min, nRows, nCols);
//
//	low  = mean - 2 * stdv;
//	high = mean + 2 * stdv;
//
//	for (int i = 0; i < nPxls; i++)
//	{
//		fPxl = (fltPixel_t)*pSrc++;
//		if (fPxl < low)
//		{
//			fPxl = low;
//		}
//		else if (fPxl > high)
//		{
//			fPxl = high;
//		}
//
//		*pBuf++ = (algPixel_t) (((fPxl - low) / (high - low)) * (float)nGrayLevels);
//	}
//
//	memcpy(pDst, pTmp, nRows * nCols * sizeof(algPixel_t));
//	free(pTmp);
//
//	return 0;
//}

int subsamp(algPixel_t *in, algPixel_t *out, int nRows, int nCols, int rowSkip, int colSkip, int *nRowsOut, int *nColsOut)
{
	algPixel_t *pSrc = in;
	algPixel_t *pDst = out;
	int row, col;
	int outRows = 0;
	int outCols = 0;

	for (row = 0; row < nRows; row += rowSkip)
	{
		outCols = 0;
		for (col = 0; col < nCols; col += colSkip)
		{
			*pDst++ = pSrc[row * nCols + col];
			outCols++;
		}
		*nColsOut = outCols;
		outRows++;
	}
	*nRowsOut = outRows;
	return 0;
}


int interp2(algPixel_t *src, algPixel_t *dst, int nSrcRows, int nSrcCols, int *nDstRows, int *nDstCols)
{
	algPixel_t *pBuf = dst;
	int row, col;

	*nDstRows = nSrcRows * 2;
	*nDstCols = nSrcCols * 2;

	for (row = 0; row < nSrcRows; row++)
	{
		for (col = 0; col < nSrcCols; col++)
		{
			*pBuf++ = src[row * nSrcCols + col];
			*pBuf++ = (src[row * nSrcCols + col] + 
				       src[row * nSrcCols + col + 1]) / 2;
		}
		for (col = 0; col < nSrcCols; col++)
		{
			*pBuf++ = (src[row * nSrcCols + col] + 
				       src[(row + 1) * nSrcCols + col]) / 2;
			*pBuf++ = (src[row * nSrcCols + col] + 
				       src[row * nSrcCols + col + 1] +
				       src[(row + 1) * nSrcCols + col] + 
					   src[(row + 1) * nSrcCols + col + 1]) / 4;
		}
	}

	return 0;
}
           

int conv2d(algPixel_t *in, algPixel_t *out, int nRows, int nCols, fltPixel_t *filter, float normFactor, int nFilterRows, int nFilterCols)
{
	float sum = 0.0;
	int m = 0, n = 0;
	int i = 0, j = 0;
	int row = 0, col = 0;
	int rowOffset = nFilterRows / 2;
	int colOffset = nFilterCols / 2;
	int rowBegIndex = rowOffset;
	int colBegIndex = colOffset;
	int pxlPos = 0;
	int fltPos = 0;

	algPixel_t *tmpBuf = (algPixel_t *)calloc((nRows + nFilterRows) * (nCols + nFilterCols), sizeof(algPixel_t));
	if (!tmpBuf)
	{
		fprintf(stderr, "File %s, Line %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	// move input to tmp buf; assumes src & dst types are the same.

	for (row = 0; row < nRows; row++)
	{
		{
			memcpy((void *)(tmpBuf + (row + rowOffset) * (nCols + nFilterCols) + colOffset), 
				   (void *)(in + row * nCols), 
				   nCols * sizeof(algPixel_t));
		}
	}
 
	for (row = rowBegIndex; row < nRows + rowOffset; row++)
	{
		for (col = colBegIndex; col < nCols + colOffset; col++)
		{
			sum = 0;
			m = 0;
			for (i = row - rowOffset; i <= row + rowOffset; i++)
			{
				n = 0;
				for (j = col - colOffset; j <= col + colOffset; j++)
				{
					pxlPos = i * (nCols + nFilterCols) + j;
					fltPos = m * nFilterCols + n;
					sum += (tmpBuf[pxlPos] * filter[fltPos]);
					n++;
				}
				m++;
			}
			out[(row - rowBegIndex) * nCols + (col - colBegIndex)] = (algPixel_t) (sum / normFactor);
		}
	}

	free((void *)tmpBuf);

	return 0;
}

// out buffer is (nrows*2) X (ncols*2)
int interpolator(algPixel_t *in, algPixel_t *out, int nRows, int nCols, fltPixel_t *filter, float normFactor, int nFilterRows, int nFilterCols)
{
	int nRows2 = nRows * 2;
	int nCols2 = nCols * 2;
	int row = 0, row2 = 0;
	int col = 0, col2 = 0;
	int err = 0;
	algPixel_t *tmpBuf = (algPixel_t *)calloc(sizeof(algPixel_t), nRows2 * nCols2);

	if (!tmpBuf)
	{
		fprintf(stderr, "File %s, Line %d, Memory Allocation Error\n", __FILE__, __LINE__);
		return -1;
	}

	// move input to tmp buf

	for (row = 0, row2 = 0; row < nRows; row++, row2 += 2)
	{
		for (col = 0, col2 = 0; col < nCols; col++, col2 += 2)
		{
			tmpBuf[row2 * nCols2 + col2] = in[row * nCols + col];
		}
	}

	err = conv2d(tmpBuf, out, nRows2, nCols2, filter, normFactor, nFilterRows, nFilterCols);

	free((void *)tmpBuf);

	return err;
}
