#!/bin/bash
#
make clean
make 

mv ./Initial ../Initial
# mv ./modelUpdate ../modelUpdate
cd ../

./Initial
#./modelUpdate
