.SUFFIXES: .tex .ps .pdf .png .jpg

all: ubercal.pdf

%.jpg: %.ps
	convert -density 150 $< $*.jpg

%.pdf: %.tex
	pdflatex $<
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
	- bash -c " ( grep Rerun $*.log && pdflatex $< ) || echo noRerun "
