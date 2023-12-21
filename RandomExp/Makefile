# Makefile

# suffixes.
.SUFFIXES: .c .o .tex .dvi .ps .pdf

# macros.
CC = gcc
CFLAGS = -O2 -Wall -lm
#CFLAGS = -O2 -Wall -lm -ggdb -lefence
DVIPS = dvips
DVIPDF = dvipdfmx -p a4
PSPDF = ps2pdf14
LATEX = platex --kanji=utf-8
PDFLATEX = pdflatex
BIBTEX = pbibtex
OBJS = null
RM = rm -f
TARGET = null

# cleaning rules.
clean:
	$(RM) *.o *~ \#*
paperclean:
	$(RM) *.o *~ \#* $(PAPER).bbl $(PAPER).blg $(PAPER).aux $(PAPER).idx $(PAPER).toc $(PAPER).lo*


# generation rules.
# PAPER = article
# BIB = article_bib
#
# $(PAPER).dvi:    $(PAPER).tex $(PAPER).bbl
#	$(BIBTEX) $(PAPER)
#	$(LATEX) $(PAPER).tex
#	$(LATEX) $(PAPER).tex
#
# $(PAPER).pdf:    $(PAPER).tex $(PAPER).bbl
#	$(BIBTEX) $(PAPER)
#	$(PDFLATEX) $(PAPER).tex
#	$(PDFLATEX) $(PAPER).tex
#
# $(PAPER).bbl:    $(BIB).bib
#	$(LATEX) $(PAPER).tex
#	$(BIBTEX) $(PAPER)



# suffix rules: "$@ is dependent on $<".
.c.o:
	$(CC) $(CFLAGS) -c $<
#.tex.dvi: $<
#	$(LATEX) $<
#	$(LATEX) $<

#.dvi.pdf: $<
#	$(DVIPDF) $<

.tex.pdf: $<
	$(PDFLATEX) $<
	$(PDFLATEX) $<


# dependencies.
