#!/bin/bash
set -e

parallel=1
mem=10000
time_=24
logfile=/cluster/work/grlab/projects/projects2020_OHSU/plots/heatmaps/collect_heatmap.log

for io in {0..542}; do 
	cmd="python heatmap_collect.py --io=${io} > ${logfile}.o 2>&1" 
	echo $cmd | bsub -n ${parallel} -J ${io}collect_heat -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -o $logfile 
done
