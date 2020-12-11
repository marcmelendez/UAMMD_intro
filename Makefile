all: tangle weave

tangle:
	txt2tangle uammd_intro.tex
	for file in `ls chapters/*.tex`; do txt2tangle $$file; done

weave:
	latex uammd_intro.tex
	dvips uammd_intro.dvi
	ps2pdf uammd_intro.ps
