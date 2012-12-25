documentation:
	doxygen Doxyfile;
	latex guide.tex; latex guide.tex; dvipdf guide.dvi;
guide:
	latex guide.tex; latex guide.tex; dvips guide.dvi; ps2pdf guide.ps

clean: 
	rm -rf *~ html latex *.aux *.bbl *.blg *.dvi *.log *.toc *.ps *.pdf;
