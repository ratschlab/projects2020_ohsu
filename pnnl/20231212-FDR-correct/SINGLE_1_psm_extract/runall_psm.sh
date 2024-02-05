#!/bin/bash

proteomics_dir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024

gitfolder=${PWD}

experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
search_res="${proteomics_dir}/OHSU/*/tide_search"
out_folder="${proteomics_dir}/OHSU"
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_psm.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder}



experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
search_res="${proteomics_dir}/ETH/*/tide_search"
out_folder="${proteomics_dir}/ETH"
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_psm.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder}
