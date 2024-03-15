#!/bin/bash

#SBATCH --job-name=TrypPepsplit
#SBATCH --output=TrypPepsplit.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=50GB

script_home=$1
list_experiments=$2
base_pipeline_folder=$3
out_folder=$4
samples=$5
make_sub_folder=$6

python ${script_home}/trypPep_to_experiments.py --list-experiments "${list_experiments}" --base-pipeline-folder "${base_pipeline_folder}" --save-folder "${out_folder}" --samples ${samples} --create-sample-subfolder ${make_sub_folder}
echo "python ${script_home}/trypPep_to_experiments.py --list-experiments "${list_experiments}" --base-pipeline-folder "${base_pipeline_folder}" --save-folder "${out_folder}" --samples ${samples} --create-sample-subfolder ${make_sub_folder}"  >> ${script_home}/run_example.sh
