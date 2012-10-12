#!/bin/sh

HomeDir=~/VFI

CPU="TRUE"
GPU=1
#MATLAB="TRUE"

if ["$GPU" -ge 1]; then
    cd $HomeDir/Thrust
    make veryclean; make; ./main
fi

#cd Matlab
#matlab -nodesktop -r main -logfile main.log
