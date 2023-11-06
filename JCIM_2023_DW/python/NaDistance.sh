#!/bin/bash

path=$(pwd)
if [ ! -d Distance ]; then
    mkdir Distance
fi
cd Distance
for first in 262;do
    for last in 274;do
	if [ ! -d dist_Na_closest_mean_${first}_${last} ]; then
	    mkdir dist_Na_closest_mean_${first}_${last}
	fi
	cd dist_Na_closest_mean_${first}_${last}

	python /home/wud18/python/SodiumLoopIonOnOffTrial.py -s ${path}/protein_Na_right_angle.pdb -t ${path}/protein_Na_stride1_aligned.dcd -sel protein\\\ and\\\ resid\\\ ${first}:${last} -title 'Closest\ mean\ distance\ between\ Na$\{^+}$\ and\ 220s\ loop' -tm 10 -nm dist_Na_closest_mean_${first}_${last} -o ${path}/Distance/dist_Na_closest_mean_${first}_${last}
	cd ..
    done
done
