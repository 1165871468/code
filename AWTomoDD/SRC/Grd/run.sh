#!/bin/bash
#
make 

mv ./Gradient ../Gradient

cd ../

mpirun -np 3 ./Gradient
