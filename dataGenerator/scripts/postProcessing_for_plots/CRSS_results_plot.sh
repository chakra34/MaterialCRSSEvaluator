#!/usr/bin/env bash

file=$1
std_size=$2                                 # precision (mean of standard deviation)
rms=${3}                                    # accuracy  (rms error between best sim and ref)
std_size=`echo "5.5/(1 + 5*${std_size})" | bc -l`         
rms=`echo "5.5/(1 + 5*${rms})" | bc -l`              # inverse scaling
gri -output ${file}.ps -p $GRI/plot \
line own \
symbol own \
symbolsize own \
colormap own \
noaxes \
clip \
tight \
aspect 2.0 \
frame 1.397 \
linear 'sys' %g 0 5 1 0 1 log 'CRSS' %g 1e-1 1e1 5 1 1 \
pos_reliability.txt 1 2 none dot ${std_size} cyan \
pos_precision.txt 1 2 none dot ${rms} magenta \
${file}_0.4.txt 1 2 solid bullet 0.15 gray \
${file}_0.3.txt 1 2 solid bullet 0.15 gray \
${file}_0.2.txt 1 2 solid bullet 0.15 gray \
${file}_0.1.txt 1 2 solid bullet 0.15 gray \
${file}_0.01.txt 1 2 solid bullet 0.15 gray \

gri_ps2pdf ${file}.ps