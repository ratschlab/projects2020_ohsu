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
run_incomplete='F'
while read f;
do
    	sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
	

	for partition in ${ms_datadir}/${sample}/*mzML*gz; do 
		part_name=$(basename ${partition}| sed 's/\.mzML\.gz//g')
		
		# Union experiments search
		searchdir=${outdir}/${sample}/${part_name}
		mkdir -p ${searchdir}
		cd ${searchdir}
		echo $searchdir
		database=${basedir}/neighbors_joint/${sample}/tide-indicies/final
#		if  ([[ "${run_incomplete}" == 'T' ]] && [[ -f ${searchdir}/*spectrumrecords.tmp ]]) || [[ "${run_incomplete}" == 'F' ]] ; then 
			echo ${searchdir}
			sbatch ${gitfolder}/script_search.sh ${crux_home} ${overwrite} ${searchdir} ${partition} ${database}
#		fi
		
		# Pipeline ETH search
		searchdir=${basedir}/ETH/${sample}/tide_search/${part_name}
		mkdir -p ${searchdir}
		cd ${searchdir}
		database=${basedir}/ETH/${sample}/neighbors/tide-indicies/final
		sbatch ${gitfolder}/script_search.sh ${crux_home} ${overwrite} ${searchdir} ${partition} ${database}


		# Pipeline OHSU search
		searchdir=${basedir}/OHSU/${sample}/tide_search/${part_name}
		mkdir -p ${searchdir}
		cd ${searchdir}
		database=${basedir}/OHSU/${sample}/neighbors/tide-indicies/final
		sbatch ${gitfolder}/script_search.sh ${crux_home} ${overwrite} ${searchdir} ${partition} ${database}
	done
done < ${fa_eth} #same fa_eth works too
