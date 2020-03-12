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


#include "imageio.h"

// Read next frame from (previously opened) file... assumes frames
// are sequential... leaves file open and file pointer positioned 
// at start of next frame in file... performs byteswap if needed

int readFrame(FILE *fp, void *image, int nPxls, int nBytesPerPxl, bool bSwap)
{
	__int64 *p64 = (__int64 *)image;
	unsigned long *p32 = (unsigned long *)image;
	unsigned short *p16 = (unsigned short *)image;
	int nPxlsRead = 0;
	int i = 0;

	if (fp == (FILE *)NULL)
	{
		fprintf(stderr, "File %s, Line %d, NULL fp passed to readFrame()\n", __FILE__, __LINE__);
		return -1;
	}

	nPxlsRead = fread(image, nBytesPerPxl, nPxls, fp);

	if (bSwap)
	{
		for (i = 0; i < nPxlsRead; i++)
		{
			if (nBytesPerPxl == sizeof(unsigned short))
			{
				*p16 = _byteswap_ushort(*p16);
				p16++;
			} else if (nBytesPerPxl == sizeof(unsigned long))
			{
				*p32 = _byteswap_ulong(*p32);
				p32++;
			}
			else if (nBytesPerPxl == sizeof(unsigned __int64))
			{
				*p64 = _byteswap_uint64(*p64);
				p64++;
			}
		}
	}

	return nPxlsRead;
}



// Supply "." for srcDir if files reside in current working directory

int readImage(void *image, char *srcDir, char *fileName, int nRows, int nCols, int nFrames, int nBytesPerPxl, bool bSwap)
{
	char *origDir = NULL;
	__int64 *p64 = (__int64 *)image;
	unsigned long *p32 = (unsigned long *)image;
	unsigned short *p16 = (unsigned short *)image;
	int nPxlsRead = 0;
	int i = 0;

	origDir = _getcwd(NULL, MAX_PATH);
	if (_chdir(srcDir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Could not change to directory=%s\n", __FILE__, __LINE__, srcDir);
		return -1;
	}

	FILE *fp = fopen(fileName, "rb");
	if (fp == (FILE *)NULL)
	{
		fprintf(stderr, "File %s, Line %d, Could not open %s for reading\n", __FILE__, __LINE__, fileName);
		return -2;
	}
	nPxlsRead = fread(image, nBytesPerPxl, nRows * nCols * nFrames, fp);
	fclose(fp);

	if (bSwap)
	{
		for (i = 0; i < nPxlsRead; i++)
		{
			if (nBytesPerPxl == sizeof(unsigned short))
			{
				*p16 = _byteswap_ushort(*p16);
				p16++;
			} else if (nBytesPerPxl == sizeof(unsigned long))
			{
				*p32 = _byteswap_ulong(*p32);
				p32++;
			}
			else if (nBytesPerPxl == sizeof(unsigned __int64))
			{
				*p64 = _byteswap_uint64(*p64);
				p64++;
			}
		}
	}

	if (_chdir(origDir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Could not change to directory=%s\n", __FILE__, __LINE__, origDir);
		return -1;
	}
	free(origDir);

	return nPxlsRead;
}



int saveFrame(void *image, char *dstDir, char *baseFileName, int nRows, int nCols, int frameNo, int nBytesPerPxl, bool bSwap)
{
	char *origDir = NULL;
	__int64 *p64 = (__int64 *)image;
	unsigned long *p32 = (unsigned long *)image;
	unsigned short *p16 = (unsigned short *)image;

	char fullFileName[MAX_PATH];
	int nPxlsToWrite = nRows * nCols;
	int i = 0;

	origDir = _getcwd(NULL, MAX_PATH);

	if (_chdir(dstDir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Could not change to directory=%s\n", __FILE__, __LINE__, dstDir);
		return -1;
	}

	sprintf(fullFileName, "%s_%dR_%dC_%dBpp_Frame%d.raw", baseFileName, nRows, nCols, nBytesPerPxl, frameNo);
	if (bSwap)
	{
		for (i = 0; i < nPxlsToWrite; i++)
		{
			if (nBytesPerPxl == sizeof(unsigned short))
			{
				*p16 = _byteswap_ushort(*p16);
				p16++;
			} else if (nBytesPerPxl == sizeof(unsigned long))
			{
				*p32 = _byteswap_ulong(*p32);
				p32++;
			}
			else if (nBytesPerPxl == sizeof(unsigned __int64))
			{
				*p64 = _byteswap_uint64(*p64);
				p64++;
			}
		}
	}
#ifdef _USE_FOPEN_FWRITES
	FILE *fp = fopen(fullFileName, "wb");
	if (fp == (FILE *)NULL)
	{
		fprintf(stderr, "File %s, Line %d, Failed fopen() on file: %s\n", __FILE__, __LINE__, fullFileName);
        return -1; 
	}
	if (fwrite((void *)image, nBytesPerPxl, nPxlsToWrite, fp) != (size_t)nPxlsToWrite)
	{
		fclose(fp);
		fprintf(stderr, "File %s, Line %d, Failed fwrite() on file: %s\n", __FILE__, __LINE__, fullFileName);
		return -1;
	}
	fclose(fp);
	if (_chdir(origDir) == -1)
	{
		fprintf(stderr, "File %s, Line %d, Could not change to directory=%s\n", __FILE__, __LINE__, origDir);
		return -1;
	}
	return nPxlsToWrite * nBytesPerPxl;
#else
	DWORD nBytesToWrite = nRows * nCols * nBytesPerPxl;
	DWORD nBytesWritten = 0;
    HANDLE hFile = CreateFileA(fullFileName,             
                        GENERIC_WRITE | GENERIC_READ,     
                        NULL, //FILE_SHARE_READ | FILE_SHARE_WRITE,	 
                        NULL,                
                        CREATE_ALWAYS,  
                        FILE_ATTRIBUTE_NORMAL,
						//FILE_FLAG_NO_BUFFERING | FILE_FLAG_WRITE_THROUGH,
					    NULL);  
 
    if (hFile == INVALID_HANDLE_VALUE) 
    { 
		char msg[256];
		sprintf(msg, "Failed CreateFileA() on file: %s, error code = %d\n", fullFileName, GetLastError());
        perror(msg);
		fflush(stdout);
        return -1; 
    }

	if (WriteFile(hFile, image, nBytesToWrite, &nBytesWritten, NULL) == 0)
    {
		char msg[256];
		sprintf(msg, "Failed WriteFile() on file: %s, error code = %d\n", fullFileName, GetLastError());
        perror(msg);
		fflush(stdout);
        return -2;
    }

	if (!FlushFileBuffers(hFile))
    {
		char msg[256];
		sprintf(msg, "Failed FlushFileBuffers() on file: %s, error code = %d\n", fullFileName, GetLastError());
        perror(msg);
		fflush(stdout);
        return -3;
    }

	if (!CloseHandle(hFile))
	{
		char msg[256];
		sprintf(msg, "Failed CloseHandle() on file: %s, error code = %d\n", fullFileName, GetLastError());
        perror(msg);
		fflush(stdout);
		return -4;
	}

	_chdir(origDir);
	free(origDir);
	return nBytesWritten;
#endif
}

