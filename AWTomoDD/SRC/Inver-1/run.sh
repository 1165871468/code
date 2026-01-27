#!/bin/bash
#
gcc -O3 Inversion.c -o Inversion -lm

mv ./Inversion ../Inversion

cd ../

./Inversion

