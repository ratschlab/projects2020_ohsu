#!/bin/bash

scripts_home=/cluster/home/prelotla/github/projects2020_immunopepper_analysis/pepQuery/tryspine_digestion/andy_lin_scripts
gitfolder=${PWD}

mode='OHSU'
#fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list.txt
fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list_fix2jx_24012024.txt

fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt
#outdir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics
outdir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics_fix2jx_24012024


if [ ${mode} == 'ETH' ]; then
    input_list=${fa_eth}
else
   input_list=${fa_ohsu}
fi

mkdir -p ${outdir}/${mode}
readme_file=${outdir}/${mode}/README
touch ${readme_file}

while read f;
do
        now=$(date)
	echo "\n Input file is $f"
	echo "Input: ${f} ${now}" >> ${readme_file}

	curF=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	
	trypsine_dir=${outdir}/${mode}/${curF}/trypsine_digest
	mkdir -p ${trypsine_dir}
	cd ${trypsine_dir}
	echo "Writing to $trypsine_dir"
	
	logdir=${trypsine_dir}/run_trypsine.log
	mkdir -p ${logdir}
	
	sbatch ${gitfolder}/script_suite.sh ${scripts_home} ${logdir} $f


done < ${input_list}
