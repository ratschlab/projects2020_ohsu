#!/bin/bash

#SBATCH --job-name=FDRsplit
#SBATCH --output=FDRsplit.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=50GB

script_home=$1
list_experiments=$2
search_res=$3
out_folder=$4
sample_pos=$5
make_sub_folder=$6

if [ -z "$make_sub_folder" ]; then
	python ${script_home}/FDR_to_experiments.py --list-experiments "${list_experiments}" --search-out-folder "${search_res}" --save-folder "${out_folder}" --sample-search-out-folder ${sample_pos}
	echo "python ${script_home}/FDR_to_experiments.py --list-experiments "${list_experiments}" --search-out-folder "${search_res}" --save-folder "${out_folder}" --sample-search-out-folder ${sample_pos}"
else
	python ${script_home}/FDR_to_experiments.py --list-experiments "${list_experiments}" --search-out-folder "${search_res}" --save-folder "${out_folder}" --sample-search-out-folder ${sample_pos} --create-sample-subfolder ${make_sub_folder}
	echo "python ${script_home}/FDR_to_experiments.py --list-experiments "${list_experiments}" --search-out-folder "${search_res}" --save-folder "${out_folder}" --sample-search-out-folder ${sample_pos} --create-sample-subfolder ${make_sub_folder}"
fi
