#!/bin/bash
#
#gcc -O3 Inversion.c -o Inversion -lm
#make clean
file=Inversion_tt_3d

make

mv ./$file ../$file

cd ../

./$file

