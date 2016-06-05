
NAME = Equations
SRCDIR = .\src

SRC = $(SRCDIR)\$(NAME).tex
TARGET = $(SRCDIR)\$(NAME).pdf

ACROBAT = "C:\Program Files (x86)\Adobe\Acrobat Reader DC\Reader\AcroRd32"

$(TARGET):

%.pdf: %.dvi
	-taskkill /IM AcroRd32.exe
	dvipdfmx $<
	start /B $(ACROBAT) $(NAME).pdf
%.dvi: %.tex
	platex  -aux-directory=$(SRCDIR) -output-directory=$(SRCDIR) $<

.PHONY: clean index
clean:
	del $(SRCDIR)\*.aux $(SRCDIR)\*.log $(SRCDIR)\*.dvi 
	del $(SRCDIR)\*.bbl $(SRCDIR)\*.blg $(SRCDIR)\*.toc

#
# process when need to update index of reference
#
index: 
	platex -output-directory=$(SRCDIR) $(SRC)

#
# process if bibliography is updated
#
bib:
	bibtex $(SRCDIR)/$(NAME)

