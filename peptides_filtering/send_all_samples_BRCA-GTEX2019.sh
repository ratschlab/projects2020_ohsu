#!/bin/bash
set -e


### Lsf and Run Parameters
mem=20000
time_=120
local_=run_cluster #"run_local"
parallel=6
#edge_or_segm=edge
suffix="commit_e5b5b51"
#suffix="commit_0a02cfe"
echo "WARNING check activation myimmuno3"

### Inputs
uniprot=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv
## Cancer Cohorts #TODO OV
base_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102
path_cancer=${base_cancer}/cohort_mutNone
sample_type=BRCA
if [ ${sample_type} == 'OV' ] ; then 
	whitelist_cancer=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_Ov_378.tsv
	target_cancer=''
	libsize_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/TCGA_Ovarian_374/cohort_mutNone/TCGA_Ovarian_374_coding_libsize75.tsv
else
	whitelist_cancer=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102.tsv
	target_cancer='' #TODO update target 
	libsize_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/TCGA_Breast_1102/cohort_mutNone/TCGA_Breast_1102_coding_libsize75.tsv
	#TODO Update gene list
        input_Junc_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_BRCA_19077
        input_annot_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/annot_BRCA_19077	
fi

## Normal Cohorts
sample_back='GTEXcore'
if [ ${sample_back} == 'GTEX' ] || [ ${sample_back} == 'GTEXcore' ] ; then 
	libsize_normal=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/GTEX_hg38_coding_libsize75.tsv
	#TODO Update gene list
	input_Segm_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Segm_GTEX2019_19077
	input_Junc_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_GTEX2019_19077

elif [ ${sample_back} == 'AllNormals' ] || [ ${sample_back} == 'matchedNormals' ]; then 
        libsize_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/GTEX2019_copy/GTEX_hg38_TCGA_All_Normals_coding_libsize75.tsv
	input_Segm_normal='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_All_Normals/cohort_mutNone_relink/Segm_TCGAnormals_19077 /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Segm_GTEX2019_19077'
        input_Junc_normal='/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_All_Normals/cohort_mutNone_relink/Junc_TCGAnormals_19077 /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_GTEX2019_19077'
fi 

if [ ${sample_back} == 'GTEX' ]; then 
	whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist
	tag_normals='Gtex'
	batch='True'
elif [ ${sample_back} == 'GTEXcore' ]; then 
	whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist_noBrain_noTestis
	tag_normals='Gtexcore'
	batch='False'
elif [ ${sample_back} == 'AllNormals' ]; then 
	#whitelist_normal=./tmp_ALL_samples
	whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX_and_TCGA_normals/sample_full_All_Normals_plus_GTEx_sample_IDs_10-2021_lib_graph_juliannelist
	tag_normals='GtexTcga'
	batch='True'
elif [ ${sample_back} == 'matchedNormals' ]; then 
	tag_normals='Matched'
	batch='False'
	if [ ${sample_type} == 'OV' ] ; then
		whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_matched_normals/GTEx_normal_samples_12-2020_Ovary_suffix.csv
	else
		whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_matched_normals/GTEX_12-2020_Breast_and_TCGA_matched_Breast.csv
	fi
fi 

## Parameters

normalizer_cancer_libsize='400000'
normalizer_normal_libsize=${normalizer_cancer_libsize}
kmer='9'

#TODO adjust parallelism 
parallelism='100'
out_partitions=1
scratch_mem=100000 # 270000 #155000
tot_batches=10

cohort_expr_lim_cancer='1'
expr_n_limit_cancer='2'
sample_expr_lim_cancer='2' #already conf2 
expr_n_limit_normal='3' 
cohort_expr_lim_normal='1'


log_dir=${base_cancer}/lsf
mkdir -p ${log_dir}
### Main 
##TODO add argument core whitelist; all normal subset and all normals with whitelist
for cohort_expr_lim_cancer in '1'; do #'0' '5'; do 
	for expr_n_limit_cancer in '1' ; do #'2' '10' 'none'; do 
		for sample_expr_lim_cancer in '2'; do # '0' ; do #other type of splicing graph
			for expr_n_limit_normal in '1' '2' '10'; do 
				for cohort_expr_lim_normal in '0' '3' '10'; do  
					while read sample; do
						for mutation_canc in ref; do 
							## Organize folders
							
							## Generate instructions
							sample_short=$(echo $sample | sed 's,\.all,,g')
							output_norm=$(dirname ${base_cancer})/filtered_backgrounds/${suffix}_${sample_back}
							output_canc=${base_cancer}/filter_${sample}/${suffix}_a_interm_cancer
							output_dir=${base_cancer}/filter_${sample}/${suffix}_${sample_back}
							output_count=${base_cancer}/filter_${sample}/${suffix}_counts
							file_count=${output_count}/G_filtered_df_${sample_short}_samp_chrt_norm_mot_unip.tsv
							mkdir -p ${output_dir}	
							mkdir -p ${output_canc}
							mkdir -p ${output_norm}
							mkdir -p ${output_count}
							logfile=${log_dir}/${sample}.cancerspec.${suffix}.log
							test_output_exist=${output_dir}/${sample}_${mutation_canc}_SampleLim${sample_expr_lim_cancer}.0CohortLim${cohort_expr_lim_cancer}.0Across${expr_n_limit_cancer}_FiltNormalsCohortlim${cohort_expr_lim_normal}.0Across${expr_n_limit_normal}.tsv						
							
							## Cmd
							if [ ! -f "${test_output_exist}/_SUCCESS" ] ; then 	
								echo $test_output_exist
								cmd="immunopepper cancerspecif --cores $parallel --mem-per-core $mem --kmer $kmer --expression-fields-c "segmentExpr" "junctionExpr" --path-cancer-matrix-edge ${input_Junc_cancer} --ids-cancer-samples "${sample}" --mut-cancer-samples ${mutation_canc} --whitelist-cancer ${whitelist_cancer} --path-cancer-libsize ${libsize_cancer} --normalizer-cancer-libsize ${normalizer_cancer_libsize} --whitelist-normal ${whitelist_normal} --path-normal-libsize ${libsize_normal} --normalizer-normal-libsize ${normalizer_normal_libsize} --output-dir $output_dir --sample-expr-support-cancer ${sample_expr_lim_cancer} --uniprot ${uniprot} --parallelism ${parallelism} --out-partitions ${out_partitions} --path-normal-matrix-segm ${input_Segm_normal} --path-normal-matrix-edge ${input_Junc_normal} --path-normal-kmer-list ${input_annot_cancer} --cohort-expr-support-norm ${cohort_expr_lim_normal} --n-samples-lim-normal ${expr_n_limit_normal} --output-count ${file_count} --interm-dir-norm ${output_norm} --interm-dir-canc ${output_canc} --tag-prefix 'G'" #TODO add back scratch for cancer? --scratch-dir 'TMPDIR'" #TODO output count remove? 
								## None case
								if [ ${expr_n_limit_cancer} != 'none' ]; then 	
									cmd1="${cmd} --cohort-expr-support-cancer ${cohort_expr_lim_cancer} --n-samples-lim-cancer ${expr_n_limit_cancer}"
								else 
									cmd1="${cmd}"
								fi 
								
								## Batch case
								if [ ${batch} == 'True' ]; then 
									cmd2="${cmd1} --tot-batches ${tot_batches} --batch-id nbtc --tag-normals ${tag_normals}nbtc"
								else
									cmd2="${cmd1} --tag-normals ${tag_normals}"
								fi
										
								cmd3="${cmd2} > ${output_dir}/${sample}.${mutation_canc}.run_cancerspecif.${suffix}.Cec${cohort_expr_lim_cancer}.Cn${expr_n_limit_cancer}.Ces.${sample_expr_lim_cancer}.Nec${cohort_expr_lim_normal}.Nn${expr_n_limit_normal}_nbtc.log 2>&1"	
								
								## Send runs matching conditions 	
								if [[ ( ${expr_n_limit_cancer} != 'none' || ${cohort_expr_lim_cancer} == '0' ) && ( ${cohort_expr_lim_normal} != '0'  ||  ${expr_n_limit_normal} == '1' ) ]] ; then 
									## Run
									if [ "$local_" = "run_local" ] ; then
										echo "running local"
										echo $cmd3
										#$cmd3
									else
										if [ ${batch} == "True" ]; then 
											for batch_id in $(seq 0 $(( $tot_batches -1))); do 
												submit=$(echo  $cmd3 | sed "s,nbtc,${batch_id},g")
												echo $submit
												echo $submit |  bsub -n ${parallel} -J ${sample_back} -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -R "rusage[scratch=$scratch_mem]" -o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
											done
										else
											echo $cmd3
											echo $cmd3 | bsub -n ${parallel} -J ${sample_back} -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -R "rusage[scratch=$scratch_mem]" -o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
										fi
									fi
								fi
							fi
						done
					done < ./tmp_samples #/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/BRCA_5samples_spladder_full.csv #./tmp_samples
				done
			done
		done
	done
done
