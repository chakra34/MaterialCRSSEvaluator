#!/usr/bin/env bash

for i in 1-1-4-1-1 1-1-4-1-2 1-1-4-2-1 1-1-4-2-2 1-2-4-1-1 1-2-4-1-2 1-2-4-2-1 1-2-4-2-2 2-1-4-1-1 2-1-4-1-2 2-1-4-2-1 2-1-4-2-2 2-2-4-1-1 2-2-4-1-2 2-2-4-2-1 2-2-4-2-2 ;
  do 
    for lambda in 0.01 0.1 0.2 0.3 0.4 ;
      do
       echo ${i};
        mkdir ${i};
       ./scripts/dataExtraction/CRSSextractor.py -o phase.${i}.rst_${lambda}.txt phase.${i}.rst --lambda slope_${i}_${lambda} --out --slope ;
       ./scripts/dataExtraction/ActivityPlots.sh ;
       ./scripts/dataExtraction/cdfPlots.sh ;
       mv ActivityPlots.pdf ${i}/ActivityPlots_${i}_${lambda}.pdf ;
       mv CumProbPlots.pdf ${i}/CumProbPlots_${i}_${lambda}.pdf ;
   done
  done