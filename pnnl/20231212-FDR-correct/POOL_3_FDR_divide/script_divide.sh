#!/bin/bash

#SBATCH --job-name=FDRdivide
#SBATCH --output=FDRdivide.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4GB

script_home=$1
file_joint=$2
map_eth_file=$3
map_ohsu_file=$4
save_folder=$5

cmd="python ${script_home}/20240220_debug_joint_proteomics_split.py --file-joint $file_joint --map-eth-file $map_eth_file --map-ohsu-file $map_ohsu_file --save-folder $save_folder"
echo ${cmd} > ${script_home}/run_example.sh
echo ${cmd}
${cmd}
