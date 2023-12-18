#!/bin/bash

#SBATCH --job-name=FDRsinglepsm
#SBATCH --output=psm.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=2GB

list_experiments=$1
search_res=$2
out_folder=$3

python psm_to_experiments.py --list-experiments ${experiment_list} --search-out-folder ${search_res} --save-folder ${out_folder} --create-sample-subfolder --rerank-psm
