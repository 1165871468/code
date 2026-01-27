#!/bin/bash
Ns=25
Nr=12

#######################
path=`pwd`
start_time=$(date +%s)

#######################
cd SRC/Model/
# make clean
make

mv Initial ../Initial
#######################
cd ../Upd/
# make clean
make
mv Update ../Update
#######################
cd ../Source/
# make clean
make
mv timeDiff ../timeDiff
mv Error ../Error
#######################
cd ../Sens
# make clean
make
mv Sensitivity ../Sensitivity
#######################
cd ../Inver
# make clean
make
mv ./Inversion ../Inversion
#######################
cd $path

ff1=Output
ff2=Output/Error
ff3=Output/Kernel
ff4=Output/Model
ff5=Output/Receivers_syn
ff6=Output/Result
ff7=Output/Snapshot_adj
ff8=Output/Snapshot_bianry_adj
ff9=Output/Snapshot_bianry_syn
ff10=Output/Snapshot_syn
ff11=Output/Source_adj

for i in $(seq 1 11)
do
    ff="ff$i"
    ff="${!ff}"
    if [ ! -d $ff ]; then
        mkdir $ff
    fi
done

for i in $(seq 1 $Ns)
do
    ff1=Output/Error/source$i
    ff2=Output/Model/Forver
    ff3=Output/Model/Inver
    ff4=Output/Receivers_syn/source$i
    ff5=Output/Snapshot_adj/source$i
    ff6=Output/Snapshot_bianry_adj/source$i
    ff7=Output/Snapshot_bianry_syn/source$i
    ff8=Output/Snapshot_syn/source$i

    for j in $(seq 2 8)
    do
        ff="ff$j"
        ff="${!ff}"
        if [ ! -d $ff ]; then
            mkdir $ff
        fi
    done
done

#######################

#######################
end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo $cost_time"s"


