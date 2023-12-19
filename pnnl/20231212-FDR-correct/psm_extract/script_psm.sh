#!/bin/bash

#SBATCH --job-name=FDRsinglepsm
#SBATCH --output=psm.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=20GB
#SBATCH --nodelist=compute-biomed-22

script_home=$1
list_experiments=$2
search_res=$3
out_folder=$4

python ${script_home}/psm_to_experiments.py --list-experiments ${experiment_list} --search-out-folder ${search_res} --save-folder ${out_folder} --create-sample-subfolder --rerank-psm
