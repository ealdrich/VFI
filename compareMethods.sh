#!/bin/sh

HomeDir=~/VFI
#HomeDir=~/Dropbox/Academics/Research/Duke/GPU/CUDA/Code/VFI/vfi

# Run baseline code
cd $HomeDir/$1
if [ $1 = Matlab ]; then
    echo
    echo =================================================
    echo Running Baseline $1 Code        
    echo =================================================
    echo
    matlab -nosplash -nodisplay -r "main;exit"
else
    echo
    echo =================================================
    echo Running Baseline $1 Code
    echo =================================================
    echo
    make veryclean; make; ./main; make veryclean
fi

# Run comparisons
for method in "${@:2}"
do
    cd $HomeDir/$method
    if [ $method = Matlab ]; then
	echo
	echo =================================================
	echo Running $method Comparison
	echo =================================================
	echo
	matlab -nosplash -nodisplay -r "main; exit"
    else
	echo
	echo =================================================
	echo Running $method Comparison
	echo =================================================
	echo
	make veryclean; make; ./main; make veryclean
    fi
    cd $HomeDir
    echo
    echo =================================================
    echo Computing Difference between $1 and $method
    echo =================================================
    echo
    #matlab -nosplash -nodisplay -r "solutionDiff('"$1"','"$method"'); exit"
    R --vanilla --slave --args $1 $method <solutionDiff.R> Rout.log
    rm Rout.log
done

echo
echo =================================================
echo Results
echo =================================================
echo
printf '%-20s %-20s %-20s %-20s\n' "Implementation" "Solution Time" "Max Abs VF Diff" "Max Abs PF Diff"
printf '%-20s %-20f %-20f %-20f\n' $1 $(cat $HomeDir/$1/solutionTime.dat) 0 0
for method in "${@:2}"
do
    printf '%-20s %-20f %-20f %-20f\n' $method $(cat $HomeDir/$method/solutionTime.dat) $(head -1 Errors_$1_$method.dat) $(tail -1 Errors_$1_$method.dat)
done


# Compare the saved results

