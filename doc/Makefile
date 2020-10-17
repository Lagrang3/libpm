flags=-interaction=nonstopmode -halt-on-error
VPATH=tex:macros:plots
PREFIX:=pm
IMAGES:=
PLOTS:=CIC-W.pdf NGP-W.pdf TSC-W.pdf PCS-W.pdf gaussian-W.pdf pk.pdf delta-W.pdf \
    wlow-gauss-cut.pdf pk-gauss-cut.pdf \
    CIC-modes-20.pdf CIC-modes-21.pdf \
    gaussian-modes-20.pdf gaussian-modes-21.pdf \
    low_pass-modes-20.pdf low_pass-modes-21.pdf \
    NGP-modes-20.pdf NGP-modes-21.pdf \
    PCS-modes-20.pdf PCS-modes-21.pdf \
    Sinc-modes-20.pdf Sinc-modes-21.pdf \
    TSC-modes-20.pdf TSC-modes-21.pdf \
    low_pass-interpolation-20.pdf \
    low_pass-interpolation-21.pdf \
    Sinc-interpolation-20.pdf \
    Sinc-interpolation-21.pdf 

PDFS:=$(PREFIX).pdf
TEXS:=intro.tex sampling.tex forces-solver.tex biblio.tex nyquist.tex force-interpolation.tex \
dynamics.tex units.tex filter-corrections.tex
MACS:=math.tex

export

default:
	$(MAKE) -C plots $(MFLAGS)
	$(MAKE) $(MFLAGS) $(PDFS)

$(PREFIX).pdf: $(PREFIX).tex $(MACS) $(IMAGES) $(TEXS) $(PLOTS) $(MACS)
	pdflatex $(flags) $<
	pdflatex $(flags) $<


clean:
	$(MAKE) -C plots $(MFLAGS) clean
	-rm $(PDFS) *.log *.aux *.out

.DELETE_ON_ERROR: $(PREFIX).pdf

.PHONY: clean default
