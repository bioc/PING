#!/bin/bash
R CMD BATCH '--args rnw="PING-PE.Rnw"' sweave.R
pdflatex PING-PE.tex
#rm *.tex
#rm *.toc
rm *.aux
rm *.log
rm *.out
rm *.Rout
