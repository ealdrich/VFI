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
	cp makefile_gpu makefile
	echo
	echo =================================================
	echo Running Baseline $1 Code
	echo =================================================
	echo
	make veryclean; make; ./main; make veryclean
	mv startTimeThrust.dat startTimeThrustGPU.dat
	mv solTimeThrust.dat solTimeThrustGPU.dat
	mv totalTimeThrust.dat totalTimeThrustGPU.dat
	mv valFunThrust.dat valFunThrustGPU.dat
	mv polFunThrust.dat polFunThrustGPU.dat
	dirName1=Thrust
    else
	if [ $1 = ThrustOMP ]; then
	    cd $HomeDir/Thrust
	    cp makefile_omp makefile
	    echo
	    echo =================================================
	    echo Running Baseline $1 Code
	    echo =================================================
	    echo
	    make veryclean; make; ./main; make veryclean
	    mv startTimeThrust.dat startTimeThrustOMP.dat
	    mv solTimeThrust.dat solTimeThrustOMP.dat
	    mv totalTimeThrust.dat totalTimeThrustOMP.dat
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
	    cp makefile_gpu makefile
	    echo
	    echo =================================================
	    echo Running $method Comparison
	    echo =================================================
	    echo
	    make veryclean; make; ./main; make veryclean
	    mv startTimeThrust.dat startTimeThrustGPU.dat
	    mv solTimeThrust.dat solTimeThrustGPU.dat
	    mv totalTimeThrust.dat totalTimeThrustGPU.dat
	    mv valFunThrust.dat valFunThrustGPU.dat
	    mv polFunThrust.dat polFunThrustGPU.dat
	    dirName2=Thrust
	else
	    if [ $method = ThrustOMP ]; then
		cd $HomeDir/Thrust
		cp makefile_omp makefile
		echo
		echo =================================================
		echo Running $method Comparison
		echo =================================================
		echo
		make veryclean; make; ./main; make veryclean
		mv startTimeThrust.dat startTimeThrustOMP.dat
		mv solTimeThrust.dat solTimeThrustOMP.dat
		mv totalTimeThrust.dat totalTimeThrustOMP.dat
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
printf '%-20s %-20s %-20s %-20s %-20s %-20s\n' "Implementation" "GPU Start Time" "Solution Time" "Total Time" "Max Abs VF Diff" "Max Abs PF Diff"
solTime=$(cat $HomeDir/$dirName1/solTime$1.dat)
if [ $1 = CUDA-C -o $1 = ThrustGPU ]; then
    startTime=$(cat $HomeDir/$dirName1/startTime$1.dat)
    totalTime=$(cat $HomeDir/$dirName1/totalTime$1.dat)
    printf '%-20s %-20f %-20f %-20f %-20f %-20f\n' $1  $startTime  $solTime $totalTime 0 0
else 
    startTime=0
    totalTime=$solTime
    printf '%-20s %-20f %-20f %-20f %-20f %-20f\n' $1  $startTime  $solTime $totalTime 0 0
fi
for method in "${@:2}"
do
    if [ $method = ThrustGPU -o $method = ThrustOMP ]; then
	dirName2=Thrust
    else
	dirName2=$method
    fi
    solTime=$(cat $HomeDir/$dirName2/solTime$method.dat)
    vfDiff=$(head -1 Errors_$1_$method.dat)
    pfDiff=$(tail -1 Errors_$1_$method.dat)
    if [ $method = CUDA-C -o $method = ThrustGPU ]; then
	startTime=$(cat $HomeDir/$dirName2/startTime$method.dat)
	totalTime=$(cat $HomeDir/$dirName2/totalTime$method.dat)
	printf '%-20s %-20f %-20f %-20f %-20f %-20f\n' $method  $startTime  $solTime $totalTime $vfDiff $pfDiff
    else
	startTime=0
	totalTime=$solTime
	printf '%-20s %-20f %-20f %-20f %-20f %-20f\n' $method  $startTime  $solTime $totalTime $vfDiff $pfDiff
    fi
done

# Clean up data files
cd $HomeDir; rm -f *.dat
cd $HomeDir/$dirName1; rm *.dat; cd $HomeDir
for method in "${@:2}"
do
    if [ $method = ThrustGPU -o $method = ThrustOMP ]; then
	dirName2=Thrust
    else
	dirName2=$method
    fi
    cd $HomeDir/$dirName2; rm -f *.dat; cd $HomeDir
done
