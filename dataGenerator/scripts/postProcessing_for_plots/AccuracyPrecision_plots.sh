#!/usr/bin/env bash

gri -output Accuracy_Random.ps -p /Users/aritrachakraborty/Documents/gri-plotting/plot_distribution.gri \
line solid \
colormap magenta \
label \
aspect 2 \
frame 1.397 \
linear 'accuracy' %g 0 1 2 5 1 \
slope_accuracy_Random.txt 2 'slope' \
cg_accuracy_Random.txt 2 'average' \
HongMei_accuracy_Random.txt 2 'max' \

gri_ps2pdf Accuracy_Random.ps 

gri -output Precision_Random.ps -p /Users/aritrachakraborty/Documents/gri-plotting/plot_distribution.gri \
line solid \
colormap cyan \
label \
noyaxis \
aspect 2 \
frame 1.397 \
linear 'precision' %g 0 1 2 5 1 \
slope_Precision_Random.txt 2 'slope' \
cg_Precision_Random.txt 2 'average' \
HongMei_Precision_Random.txt 2 'max' \

gri_ps2pdf Precision_Random.ps 
