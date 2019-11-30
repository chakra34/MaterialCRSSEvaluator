#!/usr/bin/env bash

for k in HongMei slope cg; 
  do 
    for i in 1-1-4-1-1 1-1-4-1-2 1-1-4-2-1 1-1-4-2-2 1-2-4-1-1 1-2-4-1-2 1-2-4-2-1 1-2-4-2-2 2-1-4-1-1 2-1-4-1-2 2-1-4-2-1 2-1-4-2-2 2-2-4-1-1 2-2-4-1-2 2-2-4-2-1 2-2-4-2-2 ;
       do  
        ./precision.py data_CRSS_results_${k}_${i} > ${k}_${i}_MeanStd.txt ; 
      done ; 
    done