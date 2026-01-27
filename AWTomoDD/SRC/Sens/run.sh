#!/bin/bash
#
make

#mpicc Sensitivity0.c -o Sensitivity -lm

mv ./Sensitivity ../Sensitivity

cd ../

# ./Sensitivity
mpirun -np 10 ./Sensitivity
