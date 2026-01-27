#!/bin/bash
#python3 topo_up.py

Ns=1
Nr=1

#path=../Output/normal/
path=../../Output/Kernel/
path=./
comp=Ux

R1=-0/200000/-40000/0
J1=X3i/1i
data1="$path"Kernel_vel_s"$Ns"_r"$Nr".dat
data2="$path"Smooth_s"$Ns"_r"$Nr".dat
data3="$path"Sensitivity_s"$Ns"_r"$Nr".dat
data4="$path"vel.dat
rr=Kernel_lam_s"$Ns"_r"$Nr"

gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset MAP_TICK_LENGTH_PRIMARY 1p
gmt gmtset MAP_ANNOT_OFFSET_PRIMARY 2p
gmt gmtset FONT_ANNOT_PRIMARY 4p,Helvetica MAP_FRAME_WIDTH 0.01p

gmt begin $rr png
	gmt makecpt -Cpolar -T-1.0/1.0/0.001  -Z -D
	gmt basemap -R$R1 -J$J1 -Ba -BWSen -X1.0i -Y6.0i

	awk '{ print $1,$2,$3*2.0e+3.0 }' $data4 | gmt surface -R$R1 -I200/200 -Ggrid2dvew.grd # 2-25  1,3-26
	#awk '{ print $1,$2,$3*1.0e+16.0 }' $data4 | gmt surface -R$R1 -I200/200 -Ggrid2dvew.grd
	# awk '{ print $1,$2,$3*5.0e+2.0 }' $data3 | gmt surface -R$R1 -I200/200 -Ggrid2dvew.grd  # 2-7  1,3-8

	gmt grdimage grid2dvew.grd

    # gmt colorbar -Ba
gmt end

rm gmt.* *.grd
