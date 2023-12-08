#!/bin/bash


crux_home=~/util/crux-4.1.Linux.x86_64/bin/crux
gitfolder=${PWD}

fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list.txt
fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt
ms_datadir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/ccell_rerun_2018/data/cptac

basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics
outdir=${basedir}/tide_search_joint
mkdir -p ${outdir}

overwrite='T'
while read f;
do
    	sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	
	database=${basedir}/neighbors_joint/${sample}/tide-indicies/final # finalDb.fasta # Do we need the indices there concatenated....?

	for partition in ${ms_datadir}/${sample}/*mzML*gz; do 
		part_name=$(basename ${partition}| sed 's/\.mzML\.gz//g')
		searchdir=${outdir}/${sample}/${part_name}
		mkdir -p ${searchdir}
		cd ${searchdir}
		echo $f
		echo "results in ${outdir}"
		sbatch ${gitfolder}/script_search.sh ${crux_home} ${overwrite} ${searchdir} ${partition} ${database}
	done
done < ${fa_eth} #same fa_eth works too
