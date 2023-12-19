#!/bin/bash


gitfolder=${PWD}

fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list.txt
fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt

basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics
outdir=${basedir}/tide_search_joint

while read f;
do
    	sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	# Union experiments search
        searchdir=${outdir}/${sample}
	sbatch ${gitfolder}/script_extract.sh ${searchdir}
		
	# Pipeline ETH search
	searchdir=${basedir}/ETH/${sample}/tide_search
	sbatch ${gitfolder}/script_extract.sh ${searchdir}

	# Pipeline OHSU search
	searchdir=${basedir}/OHSU/${sample}/tide_search
	sbatch ${gitfolder}/script_extract.sh ${searchdir}

done < ${fa_eth} #same fa_eth works too
