documentation:
	doxygen Doxyfile;
	pdflatex guide; bibtex guide; pdflatex guide; pdflatex guide;
	pdflatex EMtheory; bibtex EMtheory; pdflatex EMtheory; pdflatex EMtheory;
guide:guide.tex
	pdflatex guide; bibtex guide; pdflatex guide; pdflatex guide;
EMtheory:EMtheory.tex PEC_integral_equations.tex Diel_integral_equations.tex Summary_int_eq_mom_solutions.tex Near_interactions.tex MLFMA_interactions.tex
	pdflatex EMtheory; bibtex EMtheory; pdflatex EMtheory; pdflatex EMtheory;
clean:
	rm -rf *~ html latex *.aux *.bbl *.blg *.dvi *.log *.toc *.ps *.pdf;
	rm -rf *.fls *.out *.bak *.fdb_latexmk *.gz;
clean_except_pdf:
	rm -rf *~ html latex *.aux *.bbl *.blg *.dvi *.log *.toc *.ps;
	rm -rf *.fls *.out *.bak *.fdb_latexmk *.gz;
documentation_from_docker:
	if [ "${shell docker images -q pumaemdoc 2> /dev/null}" = "" ]; then \
		docker build . -t pumaemdoc; \
	fi;
	docker run --rm -u $(shell id -u ${USER} ):$(shell id -g ${USER} ) \
       -v $(shell pwd):/opt/share -w /opt/share pumaemdoc /bin/bash -c \
       "doxygen -u && make documentation"
build_docker_doc_image:
	docker build . -t pumaemdoc
.PHONY: clean clean_except_pdf