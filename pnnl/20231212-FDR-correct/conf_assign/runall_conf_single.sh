#!/bin/bash

gitfolder=${PWD}
crux_home=/cluster/home/prelotla/util/crux-4.1.Linux.x86_64/bin/crux
overwrite='T'

experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
exp_res='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/OHSU/*/tide_search_per_experiment'
out_folder='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/OHSU'

for sample_out in $(ls -d ${exp_res}); do 
	sample=$(echo ${sample_out} | cut -f9 -d '/')
	echo $sample
	outfolder="/cluster/work/grlab/projects/projects2020_OHSU/proteomics/OHSU/${sample}/assign_conf_per_experiment"
        echo ${outfolder}
	mkdir -p ${outfolder}
    for experiment in ${sample_out}/*; do
        exp_ID=$(basename $experiment | cut -f2 -d '-'| sed 's/\.txt//g')
        conf_folder=${outfolder}/${exp_ID}
	mkdir -p ${conf_folder}
	cd ${conf_folder}
        sbatch ${gitfolder}/script_confidence.sh ${crux_home} ${conf_folder} ${experiment} ${overwrite}
   done
done



#experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
#search_res='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/ETH/*/tide_search'
#out_folder='/cluster/work/grlab/projects/projects2020_OHSU/proteomics/ETH'
#cd ${out_folder}
#echo ${out_folder}
#sbatch ${gitfolder}/script_psm.sh ${gitfolder} ${experiment_list} "${search_res}" ${out_folder}
