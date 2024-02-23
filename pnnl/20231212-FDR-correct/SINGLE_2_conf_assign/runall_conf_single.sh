#!/bin/bash

gitfolder=${PWD}

### Crema-specific parameters
crema_script=/cluster/home/prelotla/github/projects2020_ohsu/pnnl/20231212-FDR-correct/FDR_scripts/script_crema.py

### Crux-specific parameters
crux_home=/cluster/home/prelotla/util/crux-4.1.Linux.x86_64/bin/crux
overwrite='T'

## General inputs
proteomics_dir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024
FDR_method='crema' # 'crux' 'crema'

for pipeline in ETH OHSU; do 
	if [[ "${pipeline}" == 'OHSU' ]]; then
		experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_experiments_per_peptides_list.txt'
	elif [[ "${pipeline}" == 'ETH' ]]; then
		experiment_list='/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_experiments_per_peptides_list.txt'
	fi
	echo "experiment file for $pipeline is: $experiment_list"

	exp_res="${proteomics_dir}/${pipeline}/*/tide_search_per_experiment"
	out_folder="${proteomics_dir}/${pipeline}"

	for sample_out in $(ls -d ${exp_res}); do 
		sample=$(echo ${sample_out} | cut -f9 -d '/')
		echo $sample
		if [[ ${FDR_method} == 'crux' ]]; then
			outfolder="${proteomics_dir}/${pipeline}/${sample}/assign_conf_per_experiment_crux"
		elif [[ ${FDR_method} == 'crema' ]]; then
		        outfolder="${proteomics_dir}/${pipeline}/${sample}/assign_conf_per_experiment_crema"
		else
		      echo "unknown ${FDR_method}"
		fi
		echo ${outfolder}
		mkdir -p ${outfolder}
	    for experiment in ${sample_out}/*; do
		exp_ID=$(basename $experiment | cut -f2 -d '-'| sed 's/\.txt//g')
		conf_folder=${outfolder}/${exp_ID}
		mkdir -p ${conf_folder}
		cd ${conf_folder}
		if [[ ${FDR_method} == 'crux' ]]; then
			sbatch ${gitfolder}/script_confidence.sh ${crux_home} ${conf_folder} ${experiment} ${overwrite}
		elif [[ ${FDR_method} == 'crema' ]]; then
		        sbatch ${gitfolder}/script_confidence_crema.sh ${crema_script} ${experiment} ${conf_folder}
		fi
	   done
	done
done
