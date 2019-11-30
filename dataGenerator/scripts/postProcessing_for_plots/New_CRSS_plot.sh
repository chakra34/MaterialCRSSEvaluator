#!/usr/bin/env bash

for k in HongMei slope cg; 
  do 
   for i in 1-1-4-1-1 1-1-4-1-2 1-1-4-2-1 1-1-4-2-2 1-2-4-1-1 1-2-4-1-2 1-2-4-2-1 1-2-4-2-2 2-1-4-1-1 2-1-4-1-2 2-1-4-2-1 2-1-4-2-2 2-2-4-1-1 2-2-4-1-2 2-2-4-2-1 2-2-4-2-2 ; 
      do  
        stdDev=`showTable -d ${k}_${i}_MeanStd.txt` ;
        ans=`awk -v stdDev="$stdDev" 'BEGIN{print stdDev * 4}'` ;
        rms=`showTable -d ${k}_${i}_RMS.txt` ;
        val=`awk -v rms="$rms" 'BEGIN{print rms * 2}'` ;
#         IFS='-' read -ra array <<< "${i}"
#         s=$(( ${array[0]} + ${array[1]} + ${array[3]} +${array[4]} ))
#         avg=`bc <<< "scale = 2; $s / 4"`
        for l in 0.01 0.1 0.2 0.3 0.4 ;
          do 
             showTable -d < CRSS_results_${k}_${i}_${l}.txt > OnlyRatio_CRSS_results_${k}_${i}_${l}.txt
            ./adding_logScaleRatio.py OnlyRatio_CRSS_results_${k}_${i}_${l}.txt ${i} > new_CRSS_results_${k}_${i}_${l}.txt
          done
        ./CRSS_results_plot.sh new_CRSS_results_${k}_${i} ${ans} ${val}
      done ; 
    done

