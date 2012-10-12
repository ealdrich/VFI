#!/bin/sh

HomeDir=~/VFI

# Choose which implementations to compare
CPU=TRUE
GPU=TRUE
MATLAB=FALSE

# Run the code in each directory
if [ $GPU = TRUE ]; then
    cd $HomeDir/Thrust
    make veryclean; make; ./main; make veryclean
fi

if [ $CPU = TRUE ]; then
    cd $HomeDir/CPP
    make veryclean; make; ./main; make veryclean
fi

if [ $MATLAB = TRUE ]; then
    cd $HomeDir/Matlab
    matlab -nodesktop -r main -logfile main.log
fi

# Compare the saved results

