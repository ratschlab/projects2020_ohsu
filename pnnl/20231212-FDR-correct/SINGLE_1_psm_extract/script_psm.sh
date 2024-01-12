#!/bin/bash

#SBATCH --job-name=FDRsinglepsm
#SBATCH --output=psm.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=72:00:00
#SBATCH --mem-per-cpu=100GB

script_home=$1
list_experiments=$2
search_res=$3
out_folder=$4

python ${script_home}/psm_to_experiments.py --list-experiments "${list_experiments}" --search-out-folder "${search_res}" --save-folder "${out_folder}" --create-sample-subfolder --rerank-psm
echo ${script_home}/psm_to_experiments.py --list-experiments "${list_experiments}" --search-out-folder "${search_res}" --save-folder "${out_folder}" --create-sample-subfolder --rerank-psm
