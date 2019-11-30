#!/usr/bin/env bash


gri -output ActivityPlots.ps -p $GRI/plot \
line own \
linewidth 1 \
symbol own \
symbolsize 0.2 \
color own \
label initial \
clip \
frame inch 1 \
log 'Schmid factor' %g 0.5 0.1 5 1 1 log 'activity' %g 1e-3 1e1 5 1 1 \
Slope_ActivationRatioDistribution_sys1.txt 1 2 ' ' solid none blue \
Slope_ActivationRatioDistribution_sys2.txt 1 2 ' ' solid none red \
Slope_ActivationRatioDistribution_sys3.txt 1 2 ' ' solid none green \
Slope_ActivationRatioDistribution_sys4.txt 1 2 ' ' solid none orange \
ActivationRatioDistribution_sys1.txt 1 2 'bas' none dot blue \
ActivationRatioDistribution_sys2.txt 1 2 'pri' none dot red \
ActivationRatioDistribution_sys3.txt 1 2 'pyr<a>' none dot green \
ActivationRatioDistribution_sys4.txt 1 2 'pyr<c+a>' none dot orange 
gri_ps2pdf ActivityPlots.ps ; 
