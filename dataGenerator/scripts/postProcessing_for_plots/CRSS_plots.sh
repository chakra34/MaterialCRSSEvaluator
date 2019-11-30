#!/usr/bin/env bash

file=$1

gri -output new_${file}.ps -p $GRI/plot \
line own \
symbol own \
color own \
noaxes \
clip \
tight \
aspect 2 \
frame 1.397 \
linear 'sys' %g 0 5 1 0 1 log 'ratio' %g 1e-1 1e1 5 1 1 \
${file}.txt 1 2 solid none LightGray \
${file}.txt 1 3 none bullet blue \
${file}.txt 1 4 none bullet red \
${file}.txt 1 5 none bullet green \
${file}.txt 1 6 none bullet orange


gri_ps2pdf new_${file}.ps