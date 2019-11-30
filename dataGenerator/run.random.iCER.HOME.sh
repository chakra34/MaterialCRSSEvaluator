#!/bin/bash -login

## wall time, nodes, ppn
#PBS -l walltime=48:00:00,nodes=1:ppn=4

### mem
#PBS -l mem=60gb

##
#PBS -j oe

### Send an email when a job is aborted, begins or ends
#PBS -m abe

###
#PBS -N cpTiM_randallY

# load necessary modules

# change to the working directory where your code is located
cd /mnt/home/zhangc43/workbench/cpTiMatrix/SlipTrace/RandomTexture_y

# set the number of OpenMP threads
export OMP_NUM_THREADS=4

# set No. of threads in DAMASK
export DAMASK_NUM_THREADS=4

# generate all possible phase configurations
# ./genPhase.py

# work through each case
# NOTE: remove phase file once done, in case we cannot finish everything in one ran
qstat -f ${PBS_JOBID}

for me in phase.*.config; do
    cp $me phase.config;
    DAMASK_spectral -g seeker.geom -l tensile.load ;
    postResults -s \
                --separation=x,y,z \
                --cr=texture,orientation,eulerangles,f,p \
                --co=resistance_slip,shearrate_slip,accumulatedshear_slip\
                --range=10 11 1 \
                seeker_tensile.spectralOut;
    addSchmidfactors --eulers eulerangles \
                     --degrees \
                     --covera 1.587 \
                     --lattice hex \
                     --force 0.0 1.0 0.0 postProc/seeker_tensile_inc100.txt ;
    cp postProc/seeker_tensile_inc100.txt ${me/.config/.rst};
    ./countSimTrace.py ${me/.config/.rst} ;
    rm $me ;
    rm seeker_tensile* ;
done
