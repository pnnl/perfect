# -*-Mode: makefile;-*-

#*BeginCopyright*************************************************************
#
#----------------------------------------------------------------------------
# Part of PERFECT Benchmark Suite (hpc.pnnl.gov/projects/PERFECT/)
#----------------------------------------------------------------------------
#
# Copyright ((c)) 2014, Battelle Memorial Institute
# Copyright ((c)) 2014, Georgia Tech Research Corporation
# All rights reserved.
#
# 1. Battelle Memorial Institute (hereinafter Battelle) and Georgia Tech
#    Research Corporation (GTRC) hereby grant permission to any person
#    or entity lawfully obtaining a copy of this software and associated
#    documentation files (hereinafter "the Software") to redistribute
#    and use the Software in source and binary forms, with or without
#    modification.  Such person or entity may use, copy, modify, merge,
#    publish, distribute, sublicense, and/or sell copies of the
#    Software, and may permit others to do so, subject to the following
#    conditions:
# 
#    * Redistributions of source code must retain the above copyright
#      notice, this list of conditions and the following disclaimers.
# 
#    * Redistributions in binary form must reproduce the above copyright
#      notice, this list of conditions and the following disclaimer in
#      the documentation and/or other materials provided with the
#      distribution.
# 
#    * Other than as used herein, neither the name Battelle Memorial
#      Institute nor Battelle may be used in any form whatsoever without
#      the express written consent of Battelle.
# 
#      Other than as used herein, neither the name Georgia Tech Research
#      Corporation nor GTRC may not be used in any form whatsoever
#      without the express written consent of GTRC.
# 
#    * Redistributions of the software in any form, and publications
#      based on work performed using the software should include the
#      following citation as a reference:
# 
#      Kevin Barker, Thomas Benson, Dan Campbell, David Ediger, Roberto
#      Gioiosa, Adolfy Hoisie, Darren Kerbyson, Joseph Manzano, Andres
#      Marquez, Leon Song, Nathan R. Tallent, and Antonino Tumeo.
#      PERFECT (Power Efficiency Revolution For Embedded Computing
#      Technologies) Benchmark Suite Manual. Pacific Northwest National
#      Laboratory and Georgia Tech Research Institute, December 2013.
#      http://hpc.pnnl.gov/projects/PERFECT/
#
# 2. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#    FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
#    BATTELLE, GTRC, OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#    INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
#    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#    SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
#    STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
#    ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
#    OF THE POSSIBILITY OF SUCH DAMAGE.
#
#*EndCopyright***************************************************************

INPUT_SIZE ?= MEDIUM

CC ?= gcc
CFLAGS ?= -Wall -g -O3

#****************************************************************************

exe1 = pa1-app-driver

libdir = .
libs   = 

$(exe1)_SRCS = instrs.c benchmark.c dwt.c dwtfusion.c ensembles.c hist.c \
    imageio.c imageops.c lapfusion.c nf.c slowmedian.c
$(exe1)_CFLAGS = -I$(libdir) -DINPUT_SIZE=INPUT_SIZE_$(INPUT_SIZE)
$(exe1)_LDFLAGS =
$(exe1)_LIBS = $(libs)
$(exe1)_LDADD = -lm

#****************************************************************************

MK_SUBDIRS =

MK_PROGRAMS_C = $(exe1)

include ../../../Makefile-template.mk

#****************************************************************************

all.local : $(libs)

check.local :
	./$(exe1)

$(libs) :
	$(MAKE) -C $(libdir)

clean.local :
	$(RM) $(exe1) ./images/output/*
