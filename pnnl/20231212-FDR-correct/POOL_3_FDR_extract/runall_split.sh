#!/bin/bash

gitfolder=${PWD}

### OHSU POOL
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
search_res='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/OHSU/*/tide_search'
out_folder='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/OHSU'
sample_pos=8
make_sub_folder='T'
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}


### ETH POOL
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
search_res='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/ETH/*/tide_search'
out_folder='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/ETH'
sample_pos=8
make_sub_folder='T'
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}


### UNION EXTRACT OHSU EXP
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
search_res='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/OHSU/*/tide_search'
out_folder='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/assign_conf_joint_to_OHSU'
sample_pos=8
make_sub_folder='F'
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}

### INUIN EXTRACT ETH EXP
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
search_res='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/ETH/*/tide_search'
out_folder='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/assign_conf_joint_to_ETH'
sample_pos=8
make_sub_folder='F'
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}



