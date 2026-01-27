#!/bin/bash
#

make clean
make

mv ./timeDiff ../timeDiff
mv ./Error ../Error

cd ../


./timeDiff
./Error
