#!/bin/bash

#rm mpi_output.txt

make clean
make

mpirun -np 35 ./gsl_example
