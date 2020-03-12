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

/**************************/
/***    UNCLASSIFIED    ***/
/**************************/


#include "benchmark.h"

// The PERFECT TAV team has left the following code almost fully intact as it was
// received by the original author mentioned below. We have added compile-time
// default settings for the INPUT_SIZE SMALL, MEDIUM, and LARGE cases to be
// consistent with the other kernels and applications, but this application can
// accept imagery of different sizes using the command line arguments described
// below. It will use the defaults listed below if no command line arguments
// are provided. The output will be written to images/output relative to the
// current working directory.

#define INPUT_SIZE_SMALL 1
#define INPUT_SIZE_MEDIUM 2
#define INPUT_SIZE_LARGE 3

#ifndef INPUT_SIZE
	#define INPUT_SIZE INPUT_SIZE_MEDIUM
#endif

static const int DEFAULT_BYTES_PER_PIXEL = 2;
static const bool DEFAULT_SWAP_BYTES = false;
// We use the same input image for both streams for simplicity. The algorithms
// are mostly image content agnostic.
#if INPUT_SIZE == INPUT_SIZE_SMALL
	static char *DEFAULT_INPUT_FILENAME1 = "../input/input_app_small.bin";
	static char *DEFAULT_INPUT_FILENAME2 = "../input/input_app_small.bin";
	static const int DEFAULT_NUM_ROWS = 480;
	static const int DEFAULT_NUM_COLS = 640;
#elif INPUT_SIZE == INPUT_SIZE_MEDIUM
	static char *DEFAULT_INPUT_FILENAME1 = "../input/input_app_medium.bin";
	static char *DEFAULT_INPUT_FILENAME2 = "../input/input_app_medium.bin";
	static const int DEFAULT_NUM_ROWS = 1080;
	static const int DEFAULT_NUM_COLS = 1920;
#elif INPUT_SIZE == INPUT_SIZE_LARGE
	static char *DEFAULT_INPUT_FILENAME1 = "../input/input_app_large.bin";
	static char *DEFAULT_INPUT_FILENAME2 = "../input/input_app_large.bin";
	static const int DEFAULT_NUM_ROWS = 2160;
	static const int DEFAULT_NUM_COLS = 3840;
#else
    #error "Unhandled value for INPUT_SIZE"
#endif


// $ benchmark sequence1.raw sequence2.raw nRows nCols nBytesPP swapBytes
//
// EX/
// $ benchmark MW_VGA.raw LW_VGA.raw 480 640 2 1
// 
// Set argv[1] .. argv[6] parameters in the project's property page
// if run from Visual Studio, otherwise specify on command line.
//
// Assumes multiple frames are stored sequentially within a file:
// - row major format
// - no gaps between images
// - no header, no footer
// - raw pixel data with specified number bytes per pixel
// 
// NOTE: 
// argv[1] = stream A file name
// argv[2] = stream B file name
// argv[3] = nRows (stream A and B must match)
// argv[4] = nCols (stream A and B must match)
// argv[5] = nBytesPerPixel (stream A and B must match)
// argv[6] = byte swap flag, 0 == no swap, 1 = swap bytes for INPUT frames
//           OUTPUT data is always written in the machine's native format
//
// For the sample images provided, byte swap is 1. 
// 
// Do NOT use this version in a WINDOWS environment. Edits will be
// required; use the previously provided version instead. The references
// to the WINDOWS variant is a carry-over from the Visual Studio version.
//
// For Linux, you can just specify the full path to the input files on the
// command line; the output will appear in the current working directory.
//
// define SAVE_IMAGES in benchmark.h and rebuild to save intermediate 
// image files for debugging purposes (don't include in benchmarks).
//
// THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY, EXPRESSED, IMPLIED, 
// OR OTHERWISE INFERRED. USE AND SUITABILITY FOR ANY PARTICULAR
// APPLICATION IS SOLELY THE RESPONSIBILITY OF THE IMPLEMENTER. 
// NO CLAIM OF SUITABILITY FOR ANY APPLICATION IS MADE.
// USE OF THIS CODE FOR ANY APPLICATION RELEASES THE AUTHOR
// AND THE US GOVT OF ANY AND ALL LIABILITY.
//
// Author: John Hodapp, USA RDECOM CERDEC NVESD
// john.hodapp@us.army.mil (NIPR), 703.704.3283 (Office)
//
// US Army RDECOM CERDEC NVESD 
// RDER-NVS-SI
// 10221 Burbeck Rd
// Ft. Belvoir, VA  22050-5806
//
// Please contact the author if there are any questions regarding the benchmarking effort
//
// NOTE: No attempt was made to optimize this code for speed. Code intended to clearly
// indicate the imaging functions required for the benchmark task. Optimization left to
// the implementer for his/her target architecture.
//



/***  GENERAL FLOW...

	Open image files A and B, allocate and initialize buffers, etc..

	While (loading frame(A) and frame(B) are successful)
	{
		Do individual benchmarks:
			ZeroMemory/memset
			"Type Convert-and-Copy" loop
			imageMat2Gray
			nf(A)
			histogram(A)
			histogramEqual(A)
			DWT(A)
			DWT_inverse(A)
			dwtFuse(A,B)
			rescaleImage(A)
			lapFuse(A,B)

		Do Ensemble processing:
			Ensemble1(A)
			Ensemble2(A,B)
			Ensemble3(A,B)
	}
	Close files A and B, clean up, exit

	
***/

int main(int argc, char *argv[])
{
	FILE *fpStreamA = NULL, *fpStreamB = NULL;
	char *streamAfn = NULL;
	char *streamBfn = NULL;
	int algErr = 0;
	int numOut = 0;
	int nRows = 0;
	int nCols = 0;
	int nBytesPP = 0;
	int frameNo = 1;
	int i = 0;
	bool swapBytes = false;
        char *current_dir = NULL, *input_dir = NULL, *output_dir = NULL;

        current_dir = calloc(1024, sizeof(*current_dir));
        input_dir = calloc(1024, sizeof(*input_dir));
        output_dir = calloc(1024, sizeof(*output_dir));

        if(current_dir == NULL || input_dir == NULL || output_dir == NULL){
           free(current_dir);
           free(input_dir);
           free(output_dir);
        }

#ifdef _WIN32
     sprintf(input_dir, "%s", IMAGE_SRC_DIR);
     sprintf(output_dir, "%s", IMAGE_DST_DIR);
#else
        if(getcwd(current_dir, 1024) == NULL){
           printf("Error in obtaining the current directory\n");
           free(current_dir); free(input_dir); free(output_dir);
           exit(11);
        }else{
	   sprintf(input_dir, "%s/%s", current_dir, IMAGE_SRC_DIR);
           sprintf(output_dir, "%s/%s", current_dir, IMAGE_DST_DIR);
        }
#endif

	printf("VC Benchmarks version: %s\n", VC_BENCH_VER);

	if (argc == 1)
	{
		streamAfn = DEFAULT_INPUT_FILENAME1;
		streamBfn = DEFAULT_INPUT_FILENAME2;
		nRows = DEFAULT_NUM_ROWS;
		nCols = DEFAULT_NUM_COLS;
		nBytesPP = DEFAULT_BYTES_PER_PIXEL;
		swapBytes = DEFAULT_SWAP_BYTES;
	}
	else if (argc == 7)
	{
		streamAfn = argv[1];
		streamBfn = argv[2];
		nRows = strtol(argv[3], NULL, 10);
		nCols = strtol(argv[4], NULL, 10);
		nBytesPP = strtol(argv[5], NULL, 10);
		swapBytes = (strtol(argv[6], NULL, 10) == 1) ? true : false;
	}
	else
	{
		fprintf(stderr, "Usage: %s [image1.raw image2.raw nRows nCols nBytesPP swapBytes]\n",
			argv[0]);
		fprintf(stderr, "\tIf no command-line arguments are given, defaults consistent with "
			"the value of INPUT_SIZE during compilation will be used.\n");
		return -1;
	}

	int nPxls = nRows * nCols;
	int nBins = (1 << (nBytesPP * 8));

	// Verify existence of destination directory....

	if (_chdir(output_dir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Destination directory (%s) error\n", __FILE__, __LINE__, output_dir);
		return -1;
	}

	if (_chdir(current_dir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Destination directory (%s) error\n", __FILE__, __LINE__, output_dir);
		return -1;
	}

#if 0
	// Switch to input directory and open up the two image files...

	if (_chdir(input_dir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Source directory (%s) error", __FILE__, __LINE__, input_dir);
		return -1;
	}
#endif

	// Open up the input files...
	
	if ((fpStreamA = fopen(streamAfn, "rb")) == (FILE *)NULL)
	{
		fprintf(stderr, "File: %s, Line: %d, Could not open %s\n", __FILE__, __LINE__, streamAfn);
		return -2;
	}
	if ((fpStreamB = fopen(streamBfn, "rb")) == (FILE *)NULL)
	{
		fclose(fpStreamA);
		fprintf(stderr, "File: %s, Line: %d, Could not open %s\n", __FILE__, __LINE__, streamBfn);
		return -2;
	}

	
	// Allocate various buffers...

	int *histA = (int *)calloc(nBins, sizeof(int));
	int *histB = (int *)calloc(nBins, sizeof(int));

	senPixel_t *rawBufA = (senPixel_t *)calloc(nPxls, sizeof(senPixel_t));
	senPixel_t *rawBufB = (senPixel_t *)calloc(nPxls, sizeof(senPixel_t));
	algPixel_t *imgSrcA = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *imgSrcB = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *buf1 = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *buf2 = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *scratch1 = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *scratch2 = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *dwtFusionOut = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *lapFusionOut = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *ensemble1Out = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *ensemble2Out = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));
	algPixel_t *ensemble3Out = (algPixel_t *)calloc(nPxls, sizeof(algPixel_t));

	if (!(rawBufA && rawBufB && imgSrcA && imgSrcB && 
		  buf1 && buf2 && histA && histB && 
		  scratch1 && scratch2 &&
		  dwtFusionOut && lapFusionOut &&
		  ensemble1Out && ensemble2Out && ensemble3Out))
	{
		free(histA);	// O/S will do this, but code checkers complain. Just call
		free(histB);	// free() on all buffers; calling free() on NULL pointer is OK.
		free(rawBufA); 
		free(rawBufB);
		free(imgSrcA); 
		free(imgSrcB);
		free(buf1); 
		free(buf2);
		free(scratch1);  
		free(scratch2);
		free(dwtFusionOut);
		free(lapFusionOut);
		free(ensemble1Out);
		free(ensemble2Out);
		free(ensemble3Out);
		fclose(fpStreamA);
		fclose(fpStreamB);
		fprintf(stderr, "File: %s, Line: %d, Memory allocation error\n", __FILE__, __LINE__);
		return -1;
	}



	// *****************************************************************************************
	// 
	// MAIN PROCESSING LOOP
	//
	// *****************************************************************************************

	while ((readFrame(fpStreamA, rawBufA, nPxls, nBytesPP, swapBytes) == nPxls)  &&
		   (readFrame(fpStreamB, rawBufB, nPxls, nBytesPP, swapBytes) == nPxls))
	{
		printf("Processing Frame # %d\n", frameNo);
		fflush(stdout);

		// Clear out buffers.... 

		memset(imgSrcA, 0, nPxls * sizeof(algPixel_t));		// <-- BENCHMARK THIS, BENCH #1
		memset(imgSrcB, 0, nPxls * sizeof(algPixel_t));
		memset(buf1, 0, nPxls * sizeof(algPixel_t));
		memset(buf2, 0, nPxls * sizeof(algPixel_t));
		memset(scratch1, 0, nPxls * sizeof(algPixel_t));
		memset(scratch2, 0, nPxls * sizeof(algPixel_t));
		memset(dwtFusionOut, 0, nPxls * sizeof(algPixel_t));
		memset(lapFusionOut, 0, nPxls * sizeof(algPixel_t));
		memset(ensemble1Out, 0, nPxls * sizeof(algPixel_t));
		memset(ensemble2Out, 0, nPxls * sizeof(algPixel_t));
		memset(ensemble3Out, 0, nPxls * sizeof(algPixel_t));
		memset(histA, 0, nBins * sizeof(int));
		memset(histB, 0, nBins * sizeof(int));


		// "Type Convert-and-Copy" raw data to buffers for processing...

		// BENCHMARK THIS LOOP: BENCH #2
		// BENCHMARK NAME: TYPE CONVERSION
		for (i = 0; i < nRows * nCols; i++)
		{
			imgSrcA[i] = (algPixel_t)rawBufA[i];
			imgSrcB[i] = (algPixel_t)rawBufB[i];
		}
		
		algErr = imageMat2Gray(imgSrcA, buf1, 16, nRows, nCols);			// <-- BENCHMARK THIS, BENCH #3 	
		assert(algErr == 0);


		// Simple noise filter (median)
		// input and output buffers can match here...
		memcpy(buf1, imgSrcA, nRows * nCols * sizeof(algPixel_t));
		memcpy(buf2, imgSrcB, nRows * nCols * sizeof(algPixel_t));

		algErr |= nf(buf1, scratch1, nRows, nCols);						// <-- BENCHMARK THIS, BENCH #4
		algErr |= nf(buf2, scratch2, nRows, nCols);					
		assert(algErr == 0);

		memcpy(buf1, scratch1, nRows * nCols * sizeof(algPixel_t));		// Put back into original buffers,
		memcpy(buf2, scratch2, nRows * nCols * sizeof(algPixel_t));		// fixes stupid error in haste.

#ifdef SAVE_IMAGES
		numOut = saveFrame(buf1, output_dir, "nf_A", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
		numOut = saveFrame(buf2, output_dir, "nf_B", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif

		// Histogram 
		algErr |= hist(imgSrcA, histA, nRows, nCols, 16);				// <-- BENCHMARK THIS, BENCH #5 
		algErr |= hist(imgSrcB, histB, nRows, nCols, 16);				
		assert(algErr == 0);

		// Histogram Equalization 
		// Assumes a histogram was previously computed...
		// 

		algErr |= histEq(imgSrcA, imgSrcA, histA, nRows, nCols, 16, 10);		// <-- BENCHMARK THIS, BENCH #6  
		algErr |= histEq(imgSrcB, imgSrcB, histB, nRows, nCols, 16, 10);	
		assert(algErr == 0);

#ifdef SAVE_IMAGES
		numOut = saveFrame(imgSrcA, output_dir, "histEq_A", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
		numOut = saveFrame(imgSrcB, output_dir, "histEq_B", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif

		// DWT Forward operation
		memcpy(buf1, imgSrcA, nPxls * sizeof(algPixel_t));
		memcpy(buf2, imgSrcB, nPxls * sizeof(algPixel_t));
		algErr |= dwt53(buf1, nRows, nCols);							// <-- BENCHMARK THIS, BENCH #7
		algErr |= dwt53(buf2, nRows, nCols);							
		assert(algErr == 0);

#ifdef SAVE_IMAGES
		numOut = saveFrame(buf1, output_dir, "dwtFwd_A", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
		numOut = saveFrame(buf2, output_dir, "dwtFwd_B", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif

		// DWT Reverse operation
		algErr |= dwt53_inverse(buf1, nRows, nCols);					// <-- BENCHMARK THIS, BENCH #8
		algErr |= dwt53_inverse(buf2, nRows, nCols);					
		assert(algErr == 0);

#ifdef SAVE_IMAGES
		numOut = saveFrame(buf1, output_dir, "dwtInv_A", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
		numOut = saveFrame(buf2, output_dir, "dwtInv_B", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif

		// DWT Fusion
		algErr |= dwtFusion(buf1, buf2, dwtFusionOut, nRows, nCols);	// <-- BENCHMARK THIS, BENCH #9
		algErr |= rescaleImage(dwtFusionOut, nRows, nCols, 16);		// <-- BENCHMARK THIS, BENCH #10 
		assert(algErr == 0);

#ifdef SAVE_IMAGES
		numOut = saveFrame(dwtFusionOut, output_dir, "dwtFused", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif

		// LAP Fusion

#ifdef _USE_DOUBLEDENSITY_LAP_FUSION	
		algErr |= lapFusionDD(buf1, buf2, lapFusionOut, nRows, nCols);	// <-- BENCHMARK THIS, BENCH #11
#else
		algErr |= lapFusion(buf1, buf2, lapFusionOut, nRows, nCols);
#endif
		assert(algErr == 0);

#ifdef SAVE_IMAGES		
		algErr |= rescaleImage(lapFusionOut, nRows, nCols, 16);	
		numOut = saveFrame(lapFusionOut, output_dir, "lapFused", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif

		// ****************
		// ENSEMBLE METRICS
		// ****************

		// Reset data...
		for (i = 0; i < nRows * nCols; i++)
		{
			imgSrcA[i] = (algPixel_t)rawBufA[i];
			imgSrcB[i] = (algPixel_t)rawBufB[i];
		}

		algErr |= ensemble_1(imgSrcA, ensemble1Out, nRows, nCols, 16, 10);			// <-- BENCHMARK THIS, BENCH #12 
		algErr |= ensemble_2(imgSrcA, imgSrcB, ensemble2Out, nRows, nCols, 16, 10);	// <-- BENCHMARK THIS, BENCH #13
		algErr |= ensemble_3(imgSrcA, imgSrcB, ensemble3Out, nRows, nCols, 16, 10);	// <-- BENCHMARK THIS, BENCH #14
		assert(algErr == 0);

#ifdef SAVE_IMAGES
		algErr |= rescaleImage(ensemble2Out, nRows, nCols, 16);		//
		algErr |= rescaleImage(ensemble3Out, nRows, nCols, 16);		//
		assert(algErr == 0);

		numOut = saveFrame(ensemble1Out, output_dir, "ensemble1", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
		numOut = saveFrame(ensemble2Out, output_dir, "ensemble2", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
		numOut = saveFrame(ensemble3Out, output_dir, "ensemble3", nRows, nCols, frameNo, sizeof(algPixel_t), false);
		assert(numOut == nRows*nCols*sizeof(algPixel_t));
#endif
 
		frameNo++;
	}

	// Clean up and exit...

	fclose(fpStreamA);
	fclose(fpStreamB);

	free(imgSrcA);
	free(imgSrcB);
	free(rawBufA);
	free(rawBufB);
	free(buf1);
	free(buf2);
	free(scratch1);
	free(scratch2);
	free(histA);
	free(histB);
	free(dwtFusionOut);
	free(lapFusionOut);
	free(ensemble1Out);
	free(ensemble2Out);
	free(ensemble3Out);

        free(current_dir);
        free(input_dir);
        free(output_dir);

	return algErr;
}


