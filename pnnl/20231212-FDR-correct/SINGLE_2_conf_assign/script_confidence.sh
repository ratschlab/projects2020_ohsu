#!/bin/bash

#SBATCH --job-name=assignConf
#SBATCH --output=conf.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=50GB

crux_home=$1
output_dir=$2
search_output=$3
overwrite=$4

echo "${crux_home} assign-confidence --output-dir ${output_dir} ${search_output} --overwrite $overwrite"
${crux_home} assign-confidence --output-dir ${output_dir} ${search_output} --overwrite $overwrite 
