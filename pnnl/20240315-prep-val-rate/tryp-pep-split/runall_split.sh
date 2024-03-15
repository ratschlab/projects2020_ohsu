#!/bin/bash


### General inputs
proteomics_dir='/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024'
gitfolder=${PWD}
rm ${gitfolder}/run_example.sh
echo "example run in ${gitfolder}/run_example.sh"
samples=$(echo 'TCGA-24-1431'  'TCGA-24-2298'  'TCGA-25-1313'  'TCGA-25-1319'  'TCGA-61-2008'  'TCGA-A2-A0D2'  'TCGA-A2-A0SX'  'TCGA-AO-A0JM'  'TCGA-BH-A18V'  'TCGA-C8-A12P')

### OHSU 
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/current_OHSU_experiments_per_peptides_list.txt' #will soon be /cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/current/OHSU_experiments_per_peptides_list.txt
base_pipeline_folder="${proteomics_dir}/OHSU"
out_folder="${proteomics_dir}/OHSU"
make_sub_folder="trypsine_digest_per_experiment"
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} ${base_pipeline_folder} ${out_folder} "${samples}" ${make_sub_folder}


### ETH 
experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list_25012024.txt'
base_pipeline_folder="${proteomics_dir}/ETH"
out_folder="${proteomics_dir}/ETH"
make_sub_folder="trypsine_digest_per_experiment"
mkdir -p ${out_folder}
cd ${out_folder}
echo ${out_folder}
sbatch ${gitfolder}/script_split.sh ${gitfolder} ${experiment_list} ${base_pipeline_folder} ${out_folder} "${samples}" ${make_sub_folder}

