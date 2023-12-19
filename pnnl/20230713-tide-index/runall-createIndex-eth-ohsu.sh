#!/bin/bash

crux_home=~/util/crux-4.1.Linux.x86_64/bin/crux
scripts_home=/cluster/home/prelotla/github/projects2020_ohsu/pnnl/andy_bin
gitfolder=${PWD}

fa_eth=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/ETH_Oct2023_data/ETH_fasta_list.txt
fa_ohsu=/cluster/work/grlab/projects/projects2020_OHSU/share_OHUS_PNLL/OHSU_Oct2023_data/OHSU_fasta_list.txt
basedir=/cluster/work/grlab/projects/projects2020_OHSU/proteomics
outdir=${basedir}/neighbors_joint
dir_irrelevant=${basedir}/index_Hsapiens
mkdir -p ${outdir}

overwrite='T'
#h_sapien_fasta=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_200707/tests/pepquery/tutorial_data/tests/data_simulated/data/pepquery_mode/uniprot-proteome_UP000005640.fasta
h_sapien_fasta=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_200707/tests/pepquery/tutorial_data/tests/data_simulated/data/pepquery_mode/hsapiens.fasta
while read f;
do
    sample=$(basename $f | cut -d '_' -f2 | cut -d '-' -f1-3 )
    ohsu_dir=${basedir}/OHSU/${sample}/trypsine_digest
    eth_dir=${basedir}/ETH/${sample}/trypsine_digest
    echo $f
	
    # Run union of OHSU and ETH pipeline
    pepsimdir=${outdir}/${sample}
    mkdir -p ${pepsimdir} 
    cd ${pepsimdir} 
    echo "results in ${pepsimdir}"
    sbatch ${gitfolder}/script_index.sh ${crux_home} ${scripts_home} ${overwrite} ${h_sapien_fasta} ${dir_irrelevant} ${eth_dir} ${ohsu_dir} 'T'
    
    # Run ETH pipeline
    pepsimdir=${basedir}/ETH/${sample}/neighbors
    mkdir -p ${pepsimdir}
    cd ${pepsimdir}
    echo "results in ${pepsimdir}"
    sbatch ${gitfolder}/script_index.sh ${crux_home} ${scripts_home} ${overwrite} ${h_sapien_fasta} ${dir_irrelevant} ${eth_dir} ${ohsu_dir} 'F'
    
    # Run OHSU pipeline
    pepsimdir=${basedir}/OHSU/${sample}/neighbors
    mkdir -p ${pepsimdir}
    cd ${pepsimdir}
    echo "results in ${pepsimdir}"
    sbatch ${gitfolder}/script_index.sh ${crux_home} ${scripts_home} ${overwrite} ${h_sapien_fasta} ${dir_irrelevant} ${ohsu_dir} ${eth_dir} 'F'
done < ${fa_eth} #same fa_eth works too
