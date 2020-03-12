# -*-Mode: makefile;-*-

#*BeginPNNLCopyright*********************************************************
#
# $HeadURL$
# $Id$
#
#***********************************************************EndPNNLCopyright*

#****************************************************************************
# $HeadURL$
#
# Nathan Tallent
#****************************************************************************

#****************************************************************************
# Input: Set the following variables
#****************************************************************************

# MK_PROGRAMS_CXX, MK_PROGRAMS_C
# MK_LIBRARIES_CXX, MK_LIBRARIES_C

# MK_SUBDIRS

# MK_TARGETS_POST: post-order dependencies w.r.t. MK_SUBDIRS 
#  all: build

MK_TARGETS_POST = \
	all \
	\
	clean distclean \
	\
	install \
	\
	check

# MK_TARGETS_PRE: pre-order dependencies w.r.t MK_SUBDIRS

MK_TARGETS_PRE = \
	info


#****************************************************************************
# Basic dependencies for targets
#****************************************************************************

# _realTargetL: target names for targets

#_realTargetSfx = .real

#_realTargetL_post = $(addsuffix $(_realTargetSfx), $(MK_TARGETS_POST))

# t <-- t.local-post <-- t.real <-- t.local-pre <-- subdirs

#----------------------------------------------------------------------------

# _localTargetL_*: target names for local targets

_localTargetSfx = .local
_localTargetSfx1  = .local1

_localTargetL_pre =  $(addsuffix $(_localTargetSfx),  $(MK_TARGETS_PRE))
_localTargetL_pre1 = $(addsuffix $(_localTargetSfx1), $(MK_TARGETS_PRE))
_localTargetL_post = $(addsuffix $(_localTargetSfx),  $(MK_TARGETS_POST))

#----------------------------------------------------------------------------

# _subdirTargetL_*: target names for all (target, subdir) pairs.

_subdirTargetSfx_post = .subdir-post
_subdirTargetSfx_pre  = .subdir-pre

_subdirTargetL_post = \
  $(foreach tgt, $(MK_TARGETS_POST), \
     $(foreach dir, $(MK_SUBDIRS), $(tgt)/$(dir)$(_subdirTargetSfx_post)))

_subdirTargetL_pre = \
  $(foreach tgt, $(MK_TARGETS_PRE), \
    $(foreach dir, $(MK_SUBDIRS), $(tgt)/$(dir)$(_subdirTargetSfx_pre)))

#----------------------------------------------------------------------------

# Post-order target definitions: A post-order target t launches its
# local target after its subdir targets.

$(MK_TARGETS_POST) : % : %$(_localTargetSfx)

$(_localTargetL_post) : %$(_localTargetSfx) : \
  $(foreach dir, $(MK_SUBDIRS), %/$(dir)$(_subdirTargetSfx_post))

$(_subdirTargetL_post) : %$(_subdirTargetSfx_post) :
ifdef DEBUG
	@echo "debug post '$@': tgt: '$(*D)' dir: '$(*F)'"
endif
	@tgt=$(*D) && dir=$(*F) \
	  && $(MAKE) -C $${dir} $${tgt}

.PHONY : $(MK_TARGETS_POST)
.PHONY : $(_localTargetL_post)
.PHONY : $(_subdirTargetL_post)

#----------------------------------------------------------------------------

# Pre-order target definitions: A pre-order target launches subdir
# targets after its local target.

$(MK_TARGETS_PRE) : % : \
  %$(_localTargetSfx1) \
  $(foreach dir, $(MK_SUBDIRS), %/$(dir)$(_subdirTargetSfx_pre))

.SECONDEXPANSION:
$(_subdirTargetL_pre) : %$(_subdirTargetSfx_pre) : $$(*D)$$(_localTargetSfx)
ifdef DEBUG
	@echo "debug pre '$@': tgt: '$(*D)' dir: '$(*F)'"
endif
	@tgt=$(*D) && dir=$(*F) \
	  && $(MAKE) -C $${dir} $${tgt}

$(_localTargetL_pre1) : %$(_localTargetSfx1) : %$(_localTargetSfx)

.PHONY : $(MK_TARGETS_PRE)
.PHONY : $(_localTargetL_pre)
.PHONY : $(_localTargetL_pre1)
.PHONY : $(_subdirTargetL_pre)


#****************************************************************************
# Compilation rules
#****************************************************************************

_sfx_cpp = .cpp
_sfx_c = .c
_sfx_cu = .cu

define _program_template_cxx
  # Note: qualify .o patterns: %.o -> %-$(1).o

  $(1)_objs = \
    $$(patsubst %$(_sfx_cpp),%-$(1).o,$$(patsubst %$(_sfx_c),%-$(1).o,$$($(1)_SRCS)))

  $(1)_objs_cpp = \
    $$(patsubst %$(_sfx_cpp),%-$(1).o,$$(filter %$(_sfx_cpp),$$($(1)_SRCS)))

  $(1)_objs_c = \
    $$(patsubst %$(_sfx_c),%-$(1).o,$$(filter %$(_sfx_c),$$($(1)_SRCS)))

  $(1) : $$($(1)_objs) $$($(1)_LIBS)
	$$(CXX) -o $$@ \
		$$($(1)_CXXFLAGS) $$(CXXFLAGS) \
		$$($(1)_LDFLAGS) $$(LDFLAGS) \
		$$^ $$($(1)_LDADD)

  $$($(1)_objs_cpp) : %-$(1).o : %$(_sfx_cpp)
	$$(CXX) -c -o $$@ $$($(1)_CXXFLAGS) $$(CXXFLAGS) $$^

  $$($(1)_objs_c) : %-$(1).o : %$(_sfx_c)
	$$(CC) -c -o $$@ $$($(1)_CFLAGS) $$(CFLAGS) $$^

  _program_objs += $$($(1)_objs)
endef

define _program_template_c
  # Note: qualify .o patterns: %.o -> %-$(1).o

  $(1)_objs = \
    $$(patsubst %$(_sfx_c),%-$(1).o,$$(filter %$(_sfx_c),$$($(1)_SRCS)))

  $(1) : $$($(1)_objs) $$($(1)_LIBS)
	$$(CC) -o $$@ \
		$$($(1)_CFLAGS) $$(CFLAGS) \
		$$($(1)_LDFLAGS) $$(LDFLAGS) \
		$$^ $$($(1)_LDADD)

  $$($(1)_objs) : %-$(1).o : %$(_sfx_c)
	$$(CC) -c -o $$@ $$($(1)_CFLAGS) $$(CFLAGS) $$^

  _program_objs += $$($(1)_objs)
endef

define _program_template_c_cu
  # Note: qualify .o patterns: %.o -> %-$(1).o

  $(1)_objs_cu = \
	$$(patsubst %$(_sfx_cu),%-$(1).o,$$(filter %$(_sfx_cu),$$($(1)_SRCS_CUDA)))

  $(1)_objs = \
	$$(patsubst %$(_sfx_c),%-$(1).o,$$(filter %$(_sfx_c),$$($(1)_SRCS)))

  $(1) : $$($(1)_objs) $$($(1)_objs_cu) $$($(1)_LIBS)
	$$(CC) -o $$@ \
		$$($(1)_CFLAGS) $$($(1)_CUDA_CFLAGS) $$(CFLAGS) \
		$$($(1)_LDFLAGS) $$(LDFLAGS) \
		$$^ $$($(1)_LDADD)

  $$($(1)_objs) : %-$(1).o : %$(_sfx_c)
	$$(CC) -c -o $$@ $$($(1)_CFLAGS) $$(CFLAGS) $$^

  $$($(1)_objs_cu) : %-$(1).o : %$(_sfx_cu)
	$$(CUDA_COMPILER) -c -o $$@ $$($(1)_CUFLAGS) $$($(1)_CUDA_CFLAGS) $$(CUFLAGS) $$^

  _program_objs += $$($(1)_objs) $$($(1)_objs_cu) 
endef

define _library_template_c_cxx
  # Note: qualify .o patterns: %.o -> %-$(1).o

  $(1)_objs = \
    $$(patsubst %$(_sfx_cpp),%-$(1).o,$$(patsubst %$(_sfx_c),%-$(1).o,$$($(1)_SRCS)))

  $(1)_objs_cpp = \
    $$(patsubst %$(_sfx_cpp),%-$(1).o,$$(filter %$(_sfx_cpp),$$($(1)_SRCS)))

  $(1)_objs_c = \
    $$(patsubst %$(_sfx_c),%-$(1).o,$$(filter %$(_sfx_c),$$($(1)_SRCS)))

  $(1) : $$($(1)_objs)
	$$(AR) rcs $$@ $$^

  $$($(1)_objs_cpp) : %-$(1).o : %$(_sfx_cpp)
	$$(CXX) -c -o $$@ $$($(1)_CXXFLAGS) $$(CXXFLAGS) $$^

  $$($(1)_objs_c) : %-$(1).o : %$(_sfx_c)
	$$(CC) -c -o $$@ $$($(1)_CFLAGS) $$(CFLAGS) $$^

  _library_objs += $$($(1)_objs)
endef


$(foreach x,$(MK_PROGRAMS_CXX),$(eval $(call _program_template_cxx,$(x))))

$(foreach x,$(MK_PROGRAMS_C),$(eval $(call _program_template_c,$(x))))

$(foreach x,$(MK_PROGRAMS_C_CU),$(eval $(call _program_template_c_cu,$(x))))

$(foreach x,$(MK_LIBRARIES_C) $(MK_LIBRARIES_CXX),$(eval $(call _library_template_c_cxx,$(x))))

#$(MK_PROGRAMS_CXX) :
#	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDADD)

#$(MK_PROGRAMS_C) :
#	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^ $(LDADD)

#%.o : %.c
#	$(CC) $(CFLAGS) -c -o $@ $^

#%.o : %.cpp
#	$(CXX) $(CXXFLAGS) -c -o $@ $^


#****************************************************************************
# Additional dependencies for targets
#****************************************************************************

all : $(MK_PROGRAMS_C) $(MK_PROGRAMS_CXX) $(MK_LIBRARIES_C) $(MK_LIBRARIES_CXX) $(MK_PROGRAMS_C_CU)

install : all
install.local : all

clean :
	for x in $(_program_objs) $(_library_objs) ; do $(RM) $${x} ; done
#	x="$(_program_objs)" && if test -n "$${x}" ; then $(RM) $${x} ; fi
#	x="$(_library_objs)" && if test -n "$${x}" ; then $(RM) $${x} ; fi

distclean : clean
	for x in $(MK_PROGRAMS_C) $(MK_PROGRAMS_C_CU) $(MK_PROGRAMS_CXX) $(MK_LIBRARIES_C) $(MK_LIBRARIES_CXX) ; do $(RM) $${x} ; done
#	x="$(MK_PROGRAMS_C)" && if test -n "$${x}" ; then $(RM) $${x} ; fi
#	x="$(MK_PROGRAMS_CXX)" && if test -n "$${x}" ; then $(RM) $${x} ; fi
#	x="$(MK_LIBRARIES_C)" && if test -n "$${x}" ; then $(RM) $${x} ; fi
#	x="$(MK_LIBRARIES_CXX)" && if test -n "$${x}" ; then $(RM) $${x} ; fi

check : all
check.local : all

help :
	@echo "$(MK_TARGETS_POST) $(MK_TARGETS_PRE)"

