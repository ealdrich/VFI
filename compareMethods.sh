#!/bin/sh

HomeDir=~/VFI
#HomeDir=~/Dropbox/Academics/Research/Duke/GPU/CUDA/Code/VFI/vfi

# Run baseline code
if [ $1 = Matlab ]; then
    cd $HomeDir/$1
    echo
    echo =================================================
    echo Running Baseline $1 Code        
    echo =================================================
    echo
    matlab -nosplash -nodisplay -r "main;exit"
    dirName1=$1
else
    if [ $1 = ThrustGPU ]; then
	cd $HomeDir/Thrust
	cp makefile_device makefile
	echo
	echo =================================================
	echo Running Baseline $1 Code
	echo =================================================
	echo
	make veryclean; make; ./main; make veryclean
	mv solTimeThrust.dat solTimeThrustGPU.dat
	mv valFunThrust.dat valFunThrustGPU.dat
	mv polFunThrust.dat polFunThrustGPU.dat
	dirName1=Thrust
    else
	if [ $1 = ThrustOMP ]; then
	    cd $HomeDir/Thrust
	    cp makefile_host makefile
	    echo
	    echo =================================================
	    echo Running Baseline $1 Code
	    echo =================================================
	    echo
	    make veryclean; make; ./main; make veryclean
	    mv solTimeThrust.dat solTimeThrustOMP.dat
	    mv valFunThrust.dat valFunThrustOMP.dat
	    mv polFunThrust.dat polFunThrustOMP.dat
	    dirName1=Thrust
	else
	    cd $HomeDir/$1
	    echo
	    echo =================================================
	    echo Running Baseline $1 Code
	    echo =================================================
	    echo
	    make veryclean; make; ./main; make veryclean
	    dirName1=$1
	fi
    fi
fi

# Run comparisons
for method in "${@:2}"
do
    if [ $method = Matlab ]; then
	cd $HomeDir/$method
	echo
	echo =================================================
	echo Running $method Comparison
	echo =================================================
	echo
	matlab -nosplash -nodisplay -r "main; exit"
	dirName2=$method
    else
	if [ $method = ThrustGPU ]; then
	    cd $HomeDir/Thrust
	    cp makefile_device makefile
	    echo
	    echo =================================================
	    echo Running $method Comparison
	    echo =================================================
	    echo
	    make veryclean; make; ./main; make veryclean
	    mv solTimeThrust.dat solTimeThrustGPU.dat
	    mv valFunThrust.dat valFunThrustGPU.dat
	    mv polFunThrust.dat polFunThrustGPU.dat
	    dirName2=Thrust
	else
	    if [ $method = ThrustOMP ]; then
		cd $HomeDir/Thrust
		cp makefile_host makefile
		echo
		echo =================================================
		echo Running $method Comparison
		echo =================================================
		echo
		make veryclean; make; ./main; make veryclean
		mv solTimeThrust.dat solTimeThrustOMP.dat
		mv valFunThrust.dat valFunThrustOMP.dat
		mv polFunThrust.dat polFunThrustOMP.dat
		dirName2=Thrust
	    else
		cd $HomeDir/$method
		echo
		echo =================================================
		echo Running $method Comparison
		echo =================================================
		echo
		make veryclean; make; ./main; make veryclean
		dirName2=$method
	    fi
	fi
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
printf '%-20s %-20f %-20f %-20f\n' $1 $(cat $HomeDir/$dirName1/solTime$1.dat) 0 0
for method in "${@:2}"
do
    if [ $method = ThrustGPU -o $method = ThrustOMP ]; then
	dirName2=Thrust
    else
	dirName2=$method
    fi
    printf '%-20s %-20f %-20f %-20f\n' $method $(cat $HomeDir/$dirName2/solTime$method.dat) $(head -1 Errors_$1_$method.dat) $(tail -1 Errors_$1_$method.dat)
done


# Compare the saved results

