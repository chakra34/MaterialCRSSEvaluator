#!/usr/bin/env bash

for k in HongMei slope cg; 
  do 
    for i in 1-1-4-1-1 1-1-4-1-2 1-1-4-2-1 1-1-4-2-2 1-2-4-1-1 1-2-4-1-2 1-2-4-2-1 1-2-4-2-2 2-1-4-1-1 2-1-4-1-2 2-1-4-2-1 2-1-4-2-2 2-2-4-1-1 2-2-4-1-2 2-2-4-2-1 2-2-4-2-2 ; 
      do  
        
        stdDev=`showTable -d ${k}_${i}_MeanStd.txt` ;
        ans=`awk -v stdDev="$stdDev" 'BEGIN{print stdDev * 4}'` ;
        rms=`showTable -d ${k}_${i}_RMS.txt` ;
        val=`awk -v rms="$rms" 'BEGIN{print rms * 2}'` ;
        ./CRSS_results_plot.sh CRSS_results_${k}_${i} ${ans} ${val}
      done ; 
    done