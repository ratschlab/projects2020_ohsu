#!/bin/bash
set -e

time_=24
mem=20000
jobname=mhc

input_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102
sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/BRCA_5samples_spladder_full_split.csv
filter_id=commit_d4aee54_GTEXcore
filter_param=_SampleLim2.0CohortLim0.0Across1_FiltNormalsGtexcoreCohortlim0.0Across0_FiltUniprot.tsv

mhc_software_path=/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/sw/netMHC-4.0
alleles_file=/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/data/Shukla_Wu_Getz_Polysolver_HLA_Types_2015.tsv

head -1 ${sample_file} > ./tmp_sample_file
logdir=./logs_${jobname}
mkdir -p ${logdir}

while read sample; do
        #allele_list=$(grep $(echo ${sample} | cut -f1-3 -d '-') ${alleles_file} |  cut -f2-31 -d '	' | sed s',	,\,HLA\-',g |  sed -E s',^,HLA\-',g)	
	sample_short=$(echo $sample | cut -f1-3 -d '-')
	allele_sub_file=/cluster/work/grlab/projects/TCGA/PanCanAtlas/peptides_neoantigen/analysis_pancan/ccell_rerun_2018/output/mhc/donor_alleles.netmhc/${sample_short}.netmhc_allele_string.txt
	allele_list=$(cat $allele_sub_file)
	echo $sample
        echo $allele_list
	for mutation_type in ref; do 
		for analysis_folder in ${input_dir}/filter_${sample}/${filter_id}/G_${sample}_${mutation_type}${filter_param};  do 
			analysis_output=$(dirname $analysis_folder)/$( basename $analysis_folder |  sed 's/\.tsv//g' )_FiltMHC.tsv
			cmd="immunopepper mhcbind --partitioned-tsv ${analysis_folder} --output-dir ${analysis_folder}  --mhc-software-path ${mhc_software_path} --bind-score-method percentile_rank --argstring '--mhc-predictor netmhc4 --input-peptides-file ${analysis_folder}/kmer_only.tsv --mhc-alleles ${allele_list} --output-csv ${analysis_output}' > ${input_dir}/${sample}_${mutation_type}_run_mhcbind.log 2>&1"
			echo $cmd
			echo ${cmd} | bsub -J ${jobname} -n 1 -R "rusage[mem=${mem}]" -W ${time_}:00 -o ./${logdir}/MHC_${filter_param} 	

		done
	done
done < ${sample_file}





