#!/bin/bash

### General inputs
FDR_method='crema' #'crux' # 'crux' 'crema'
proteomics_dir='/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024'
gitfolder=${PWD}
rm ${gitfolder}/run_example.sh
echo "example run in ${gitfolder}/run_example.sh"

fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list_fixMerge_25012024.txt
fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt


### Crema thresholds

if [[ ${FDR_method} == 'crux' ]]; then
	save_suffix='FDRcrux'
	input_suffix="FDRcrux/assign-confidence.target.txt"
elif [[ ${FDR_method} == 'crema' ]]; then
	save_suffix="FDRcrema"
	input_suffix="FDRcrema/crema.peptides.txt"
else
	echo "unknown ${FDR_method}"
fi



while read f;
do
        sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	### UNION SPLIT OHSU EXP
	file_joint="${proteomics_dir}/tide_search_joint/${sample}/${input_suffix}"
	map_ohsu_file="${proteomics_dir}/OHSU/${sample}/trypsine_digest/pepID_joint_original.tsv.gz"
	map_eth_file="${proteomics_dir}/ETH/${sample}/trypsine_digest/pepID_joint_original.tsv.gz"
	save_folder="${proteomics_dir}/tide_search_joint/${sample}/${save_suffix}"

	cd ${save_folder}
	echo ${out_folder}
	sbatch ${gitfolder}/script_divide.sh ${gitfolder} ${file_joint} ${map_eth_file} ${map_ohsu_file} ${save_folder} 

done < ${fa_eth}
