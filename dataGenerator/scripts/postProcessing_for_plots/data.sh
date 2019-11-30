#!/usr/bin/env bash

for k in HongMei slope cg; 
 do 
  for i in 1-1-4-1-1 1-1-4-1-2 1-1-4-2-1 1-1-4-2-2 1-2-4-1-1 1-2-4-1-2 1-2-4-2-1 1-2-4-2-2 2-1-4-1-1 2-1-4-1-2 2-1-4-2-1 2-1-4-2-2 2-2-4-1-1 2-2-4-1-2 2-2-4-2-1 2-2-4-2-2 ; 
   do 
    for j in 0.01 0.1 0.2 0.3 0.4; 
     do 
      showTable -d CRSS_results_${k}_${i}_${j}.txt > data_CRSS_results_${k}_${i}_${j}.txt ; 
     done ; 
    done ; 
  done