flags=-interaction=nonstopmode -halt-on-error
VPATH=tex:macros:plots
PREFIX:=pm
IMAGES:=
PLOTS:=CIC-W.pdf NGP-W.pdf TSC-W.pdf PCS-W.pdf gaussian-W.pdf pk.pdf delta-W.pdf \
    wlow-gauss-cut.pdf pk-gauss-cut.pdf
PDFS:=$(PREFIX).pdf
TEXS:=intro.tex sampling.tex forces-solver.tex biblio.tex nyquist.tex force-interpolation.tex \
dynamics.tex units.tex
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
