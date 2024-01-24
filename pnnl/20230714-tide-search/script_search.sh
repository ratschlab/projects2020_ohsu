#!/bin/bash

#SBATCH --job-name=cruxsearch
#SBATCH --output=search.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=50GB
#SBATCH --exclude=compute-biomed-10
crux_home=$1
overwrite=$2
output_dir=$3
ms_partition=$4
database=$5

bashCommand="${crux_home} tide-search --overwrite T --precursor-window 40 --concat T --top-match 1000000000 --num-threads 1 --output-dir ${output_dir} ${ms_partition} ${database}"
$bashCommand
