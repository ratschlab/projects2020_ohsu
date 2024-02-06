#!/bin/bash

### General inputs
FDR_method='crux' # 'crux' 'crema'
proteomics_dir='/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024'
gitfolder=${PWD}
rm ${gitfolder}/run_example.sh
echo "example run in ${gitfolder}/run_example.sh"

### Crema thresholds

if [[ ${FDR_method} == 'crux' ]]; then
	save_suffix='_crux'
	input_suffix="FDRcrux/assign-confidence.target.txt"
elif [[ ${FDR_method} == 'crema' ]]; then
	save_suffix="_crema"
	input_suffix="FDRcrema/crema.peptides.txt"
else
	echo "unknown ${FDR_method}"
fi

#### OHSU POOL
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
search_res="${proteomics_dir}/OHSU/*/tide_search/${input_suffix}"
out_folder="${proteomics_dir}/OHSU"
sample_pos=8
make_sub_folder="assign_conf_pooled_FDR${save_suffix}"
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}


### ETH POOL
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
search_res="${proteomics_dir}/ETH/*/tide_search/${input_suffix}"
out_folder="${proteomics_dir}/ETH"
sample_pos=8
make_sub_folder="assign_conf_pooled_FDR${save_suffix}"
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}


### UNION EXTRACT OHSU EXP
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
search_res="${proteomics_dir}/tide_search_joint/*/${input_suffix}"
out_folder="${proteomics_dir}/assign_conf_joint_to_OHSU${save_suffix}"
sample_pos=8
make_sub_folder=''
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}

### INUIN EXTRACT ETH EXP
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
search_res="${proteomics_dir}/tide_search_joint/*/${input_suffix}"
out_folder="${proteomics_dir}/assign_conf_joint_to_ETH${save_suffix}"
sample_pos=8
make_sub_folder=''
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder} ${sample_pos} ${make_sub_folder}



