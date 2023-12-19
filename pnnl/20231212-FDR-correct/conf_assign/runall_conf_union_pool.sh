#!/bin/bash


gitfolder=${PWD}
crux_home=/cluster/home/prelotla/util/crux-4.1.Linux.x86_64/bin/crux
fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list.txt
fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt

basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics
outdir=${basedir}/tide_search_joint

while read f;
do
    	sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	# Union experiments search
        searchdir=${outdir}/${sample}
	search_output=${searchdir}/tide-search-concat.txt
	cd $searchdir
	sbatch script_confidence.sh ${crux_home} ${output_dir} ${search_output}
		
	# Pipeline ETH search (Pool)
	searchdir=${basedir}/ETH/${sample}/tide_search
	search_output=${searchdir}/tide-search-concat.txt
	cd $searchdir
	sbatch script_confidence.sh ${crux_home} ${output_dir} ${search_output}

	# Pipeline OHSU search (Pool)
	searchdir=${basedir}/OHSU/${sample}/tide_search
	search_output=${searchdir}/tide-search-concat.txt
	cd $searchdir
	sbatch script_confidence.sh  ${crux_home} ${output_dir} ${search_output}
done < ${fa_eth} #same fa_eth works too
~
