#!/bin/bash

#path=../Output/normal/
path=./

comp=Ux

R1=0/200000/-40000/0
J1=X5i/1i
rgb=rgb_gji.cpt
rgb=cpt/PiYG_11.cpt
rgb=seis
rgb=tong2.cpt
data="$path"vel_initial.dat
data="$path"vel.dat
# data=../Result/vel_per_in"$num".dat

rr=vel

gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_TICK_LENGTH_PRIMARY 1p
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 2p
gmt gmtset FONT_ANNOT_PRIMARY 4p,Helvetica MAP_FRAME_WIDTH 0.01p

#gmt makecpt -C$rgb -T-1.0/1.0/0.1 -D > test.cpt

gmt begin $rr png
	gmt makecpt -C$rgb -T6.0/7.5 -D

	gmt basemap -R$R1 -J$J1 -Ba -BWSen -X1.0i -Y6.0i

	#awk '{ print $1, -$2, $3/1000 }' $data | gmt surface -R$R1 -I200/200 -Ggrid2dvew.grd
	awk '{ print $1, -$2, $3*1.0e+22 }' $data | gmt surface -R$R1 -I200/200 -Ggrid2dvew.grd

	gmt grdimage grid2dvew.grd

    gmt colorbar -Ba
gmt end

rm gmt.* *.grd
