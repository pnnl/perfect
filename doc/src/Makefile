# $Id$

#****************************************************************************

# TEX setup:
#   - end with a ':' to append standard search path
#   - double trailing slash '//' searches the directory recursively
export TEXINPUTS = .:./texmf/tex/latex//:
export BSTINPUTS = .:./texmf/bibtex/bst//:
export BIBINPUTS = .:./texmf/bibtex/bib:

#****************************************************************************

TARGETS = PERFECT-Suite-manual

TARGETS_PDF = $(addsuffix .pdf,  $(TARGETS))

TARGETS_DVI = $(addsuffix .dvi,  $(TARGETS))

#****************************************************************************

# You want latexmk to *always* run, because make does not have all the info.
# Also, include non-file targets in .PHONY so they are run regardless of any
# file of the given name existing.
.PHONY: $(TARGETS_PDF) $(TARGETS_DVI) all clean

# The first rule in a Makefile is the one executed by default ("make"). It
# should always be the "all" rule, so that "make" and "make all" are identical.
all: $(TARGETS_PDF)

# CUSTOM BUILD RULES

# In case you didn't know, '$@' is a variable holding the name of the target,
# and '$<' is a variable holding the (first) dependency of a rule.
# "raw2tex" and "dat2tex" are just placeholders for whatever custom steps
# you might have.

#%.tex: %.raw
#        ./raw2tex $< > $@

#%.tex: %.dat
#        ./dat2tex $< > $@

# MAIN LATEXMK RULE

# -pdf tells latexmk to generate PDF directly (instead of DVI).
# -pdflatex="" tells latexmk to call a specific backend with specific options.
# -use-make tells latexmk to call make for generating missing files.

# -interactive=nonstopmode keeps the pdflatex backend from stopping at a
# missing file reference and interactively asking you for an alternative.

$(TARGETS_PDF) : %.pdf : %.tex
	latexmk -pdf -pdflatex="pdflatex -interactive=nonstopmode" \
	  -jobname=$* \
	  $<

$(TARGETS_DVI) : %.dvi : %.tex
	latexmk -dvi \
	  -use-make \
	  -jobname=$* \
	  $<

clean :
	latexmk -C -pdf -bibtex
