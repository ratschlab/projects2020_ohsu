#!/bin/bash

gitfolder=${PWD}

### Crema-specific parameters
crema_script=/cluster/home/prelotla/github/projects2020_ohsu/pnnl/20231212-FDR-correct/FDR_scripts/script_crema.py
eval_fdr=0.05
threshold=0.05

### Crux-specific parameters
crux_home=/cluster/home/prelotla/util/crux-4.1.Linux.x86_64/bin/crux
overwrite='T'

### General Inputs
fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list.txt
#fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list_fix2jx_24012024.txt
fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list_fixMerge_25012024.txt
fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt

#basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics
#basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fix2jx_24012024
basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fixMerge_25012024
outdir=${basedir}/tide_search_joint
FDR_method='crema' # 'crux' 'crema'
while read f;
do
    	sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	# Union experiments search
        searchdir=${outdir}/${sample}
	search_output=${searchdir}/tide-search-concat.txt
	cd $searchdir
	if [[ ${FDR_method} == 'crux' ]]; then 
		sbatch ${gitfolder}/script_confidence.sh ${crux_home} ${searchdir} ${search_output} ${overwrite}
	elif [[ ${FDR_method} == 'crema' ]]; then
	        FDR_dir=${searchdir}/FDRcrema_${threshold}
		mkdir -p ${FDR_dir}
		sbatch ${gitfolder}/script_confidence_crema.sh ${crema_script} ${search_output} ${eval_fdr} ${threshold} ${FDR_dir} 
	else
		echo "unknown ${FDR_method}"
	fi

	# Pipeline ETH search (Pool)
	searchdir=${basedir}/ETH/${sample}/tide_search
	search_output=${searchdir}/tide-search-concat.txt
	cd $searchdir
	if [[ ${FDR_method} == 'crux' ]]; then
		sbatch ${gitfolder}/script_confidence.sh ${crux_home} ${searchdir} ${search_output} ${overwrite}
	elif [[ ${FDR_method} == 'crema' ]]; then
	        FDR_dir=${searchdir}/FDRcrema_${threshold}
		mkdir -p ${FDR_dir}
		sbatch ${gitfolder}/script_confidence_crema.sh ${crema_script} ${search_output} ${eval_fdr} ${threshold} ${FDR_dir}
	else
		echo "unknown ${FDR_method}"
	fi

	# Pipeline OHSU search (Pool)
	searchdir=${basedir}/OHSU/${sample}/tide_search
	search_output=${searchdir}/tide-search-concat.txt
	cd $searchdir
	if [[ ${FDR_method} == 'crux' ]]; then
		sbatch ${gitfolder}/script_confidence.sh ${crux_home} ${searchdir} ${search_output} ${overwrite}
	elif [[ ${FDR_method} == 'crema' ]]; then
		FDR_dir=${searchdir}/FDRcrema_${threshold}
		mkdir -p ${FDR_dir}
		sbatch ${gitfolder}/script_confidence_crema.sh ${crema_script} ${search_output} ${eval_fdr} ${threshold} ${FDR_dir}
	else
		echo "unknown ${FDR_method}"
	fi
done < ${fa_eth} #same fa_eth works too
