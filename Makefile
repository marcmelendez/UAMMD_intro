all: tangle weave

tangle:
	txt2tangle uammd_intro.tex
	for file in `ls chapters/*.tex`; do txt2tangle $$file; done

weave:
	latex uammd_intro.tex
	dvips -o uammd_intro.tmp.ps uammd_intro.dvi || exit 1
	sed '/^SDict begin \[$$/ , /^end$$/d' uammd_intro.tmp.ps > uammd_intro.ps
	rm uammd_intro.tmp.ps
	ps2pdf uammd_intro.ps
