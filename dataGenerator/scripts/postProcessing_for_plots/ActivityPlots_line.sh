#!/usr/bin/env bash


gri -output ActivityPlots_all_cpTi.ps -p $GRI/plot \
line own \
linewidth 1 \
symbol own \
symbolsize 0.2 \
color own \
label initial \
clip \
frame inch 1 \
log 'SF' %g 0.5 0.1 5 1 1 log 'Activity' %g 1e-3 1e1 5 1 1 \
Slope_ActivationRatioDistribution_sys1.txt 1 2 ' ' solid none blue \
Slope_ActivationRatioDistribution_sys2.txt 1 2 ' ' solid none red \
Slope_ActivationRatioDistribution_sys3.txt 1 2 ' ' solid none green \
Slope_ActivationRatioDistribution_sys4.txt 1 2 ' ' solid none orange \
ActivationRatioDistribution_sys1.txt 1 2 '' none dot blue \
ActivationRatioDistribution_sys2.txt 1 2 '' none dot red \
ActivationRatioDistribution_sys3.txt 1 2 '' none dot green \
ActivationRatioDistribution_sys4.txt 1 2 '' none dot orange 

gri_ps2pdf ActivityPlots_all_cpTi.ps ; open ActivityPlots_all_cpTi.pdf
