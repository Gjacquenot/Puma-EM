documentation:
	doxygen Doxyfile;
	latex guide; bibtex guide; latex guide; latex guide; dvipdf guide.dvi
	latex EMtheory; bibtex EMtheory; latex EMtheory; latex EMtheory; dvipdf EMtheory.dvi
guide:
	latex guide.tex; bibtex guide; latex guide; latex guide.tex; dvips guide.dvi; ps2pdf guide.ps
EMtheory:
	latex EMtheory; bibtex EMtheory; latex EMtheory; latex EMtheory; dvips EMtheory.dvi; ps2pdf EMtheory.ps
clean: 
	rm -rf *~ html latex *.aux *.bbl *.blg *.dvi *.log *.toc *.ps *.pdf;
