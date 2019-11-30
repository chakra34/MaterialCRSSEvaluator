#!/usr/bin/env bash

obs="observation"
theo="theoretical"
#-------------- sys 1 - basal; 2 - prism; 3 - pyr<a>; 4 - pyr<c+a> -------------------- #
gri -output CumProbPlots.ps -p $GRI/plot \
line solid \
symbol own \
color own \
clip frame 5 \
log 'SF' %g 0.5 0.1 5 1 1 linear 'cumProb' %g 0 1.0 5 1 1 \
SF_${obs}Distribution_sys1.txt 1 2 dot blue \
SF_${theo}Distribution_sys1.txt 1 2 none blue \
SF_${obs}Distribution_sys2.txt 1 2 dot red \
SF_${theo}Distribution_sys2.txt 1 2 none red \
SF_${obs}Distribution_sys3.txt 1 2 dot green \
SF_${theo}Distribution_sys3.txt 1 2 none green \
SF_${obs}Distribution_sys4.txt 1 2 dot yellow \
SF_${theo}Distribution_sys4.txt 1 2 none yellow 
gri_ps2pdf CumProbPlots.ps ; 

