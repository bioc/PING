#!/bin/bash
R CMD BATCH '--args rnw="PING.Rnw"' sweave.R
pdflatex PING.tex
#rm *.tex
#rm *.toc
rm *.aux
rm *.log
rm *.out
