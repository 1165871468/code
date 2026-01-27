#!/bin/bash
num=40

#######################
path=`pwd`

start_time=$(date +%s)
#######################
cd $path"/SRC"

echo Initialize...
./Initial

echo Calculating forward wavefiled ...
mpirun -np $num ./Update

echo Sorting error ...
./timeDiff
./Error

echo Calculating sensitivity matrix ...
mpirun -np $num ./Sensitivity

echo Inversing...
./Inversion

# cd ../
######################
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo $cost_time"s"
######################
