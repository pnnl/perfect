# -*-Mode: makefile;-*-

#*BeginCopyright*************************************************************
#
# $HeadURL$
# $Id$
#
#***************************************************************EndCopyright*

#****************************************************************************
# $HeadURL$
#****************************************************************************

MK_SUBDIRS = \
	suite

include Makefile-template.mk

#****************************************************************************

TAR = $(if $(findstring Darwin, $(shell uname)), tar, tar) # gnutar
TARNM = perfect-suite-v1.0

dist :
	nm_cur=`basename ${PWD}` ; \
	nm_new=$(TARNM) ; \
	cd .. ; \
	if [[ ! -e $${nm_new} ]] ; then ln -s $${nm_cur} $${nm_new} ; fi ; \
	${TAR} zcvf $${nm_new}.tar.gz \
	  --dereference --exclude=".svn" \
	  --exclude="suite/sar/tools" \
	  --exclude="doc/src" \
	  $${nm_new}
