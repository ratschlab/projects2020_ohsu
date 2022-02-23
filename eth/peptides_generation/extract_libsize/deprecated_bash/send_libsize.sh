#!/bin/bash
set -e

time_=24
mem=20000
for folder in TCGA_Ovarian_374 TCGA_Breast_1102 TCGA_All_Normals; do 

	gene_expre_file=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/${folder}/cohort_mutNone/gene_expression_detail.pq
	result_file=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/${folder}/cohort_mutNone/${folder}_library_size_75perce.tsv


	echo "sh compute_libsizes.sh ${gene_expre_file} ${result_file} "| bsub -J libsizes -W ${time_}:00 -R "rusage[mem=${mem}]" -o /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1/${folder}/cohort_mutNone/libsize.lsf 

done
