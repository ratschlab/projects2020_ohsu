#!/bin/bash
set -e


### Lsf and Run Parameters: Can be changed as wished to submit the job to the cluster. Basic running stage 
mem=20000
time_=04
#local_="print_only_the_command" #run_cluster # Use print_only_the_command to only print the command and preview what will happen, and use run_cluster to submit to the lsf system via bsub
local_="run_cluster"
parallel=2
suffix="run_test1_commit_bb" #Choose your run name (output folder name)
echo "WARNING check activation myimmuno3"
base_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/learning_filter
mkdir -p ${base_cancer}
log_dir=${base_cancer}/lsf
mkdir -p ${log_dir}

### Inputs

## Entire proteome 
uniprot=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv

## Cancer Cohorts 
# There are two possible cancer cohorts BRCA and OV. Here we will work with BRCA cohort 
# whitelist
whitelist_cancer=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102.tsv
# libsize
libsize_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_libsizes_conf2_annotFrame_cap1000_runs_pya0.17.1_KEEP/TCGA_Breast_1102/cohort_mutNone/TCGA_Breast_1102_coding_libsize75.tsv
# Here we applied a "trick". Instead of using the real paths of the files we did use a simililink to access the files that we care about in the same place. 
input_Junc_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_BRCA_toy_data #_19077
input_annot_cancer=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/annot_BRCA_toy_data #_19077	

## Normal Cohorts
sample_back='GTEXcore'
# We work with 4 different cohorts All Normals (TCGA Normals + GTEX), GTEX, GTEXcore (same as GTEX but excluding brain and testis", Matched Normals (Only the normals for Breast tissue in the case of Breast cancer". Here we are going to focus on GTEXcore. 
# libsize
libsize_normal=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_librarysize_pya.0.17.1_conf2_annot_ref_chrall_cap1000/cohort_mutNone/GTEX_hg38_coding_libsize75.tsv
# Here we applied the same "trick" as above. We have relinked the files of interest with simililinks
input_Segm_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Segm_GTEX2019_toy_data #_19077
input_Junc_normal=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102/cohort_mutNone_relink/Junc_GTEX2019_toy_data #_19077
# whitelist
whitelist_normal=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist_noBrain_noTestis
# Name tag
tag_normals='Gtexcore'
# For very big cohorts we could use batches
batch='False'

## Normalizers
normalizer_cancer_libsize='400000'
normalizer_normal_libsize=${normalizer_cancer_libsize}
kmer='9'

## Spark Specific parameters: Can be tuned to make it quicker. But discuss first. Expert runninng stage 
parallelism='1000'
out_partitions=1
scratch_mem=2700  #270000  #270000 #155000 # 270000 #155000
tot_batches=1

### Filtering specific parameters 
# In the real code the parameters are specified in a loop. Here we will skip the loop for learning purposes. The loop was a convenient way to run many parameters in combination 
# The filtering uses a few parameters that one can set. Basic running stage. 
# First it is important to understand what they mean
expr_nsamples_limit_normal=3,10 #'Any,2' '10,2, '0,0'
# These are the two numbers passed to --cohort-expr-support-norm, --n-samples-lim-normal. I just specify the two numbers here because in the any case, one of the two arguments is ommited (skipped). 
# Above we have the "tolerance" to include a kmer in the normal background set. The first number is the number of reads (expression) and the second is the number of samples. They are specified together, however they are applied sequencially.
# Take 3, 10: This means that you are going to filter with (a) The kmer needs 3 reads in any of the normal samples. (b) The kmer needs 10 samples with any number of reads. We can re-write (a) as 3 reads in >=1 samples (b) >0 reads in >=10 samples. 
# Take Any, 2:  This means that you are going to filter with (a) -> basically nothing happens there. (b) The kmer needs 2 samples with any number of reads. 
# Now the special case 0, 0:  This means that you are going to filter with (a) The kmer needs >0 reads in any of the normal samples. (b) The kmer needs more than 0 samples with any number of reads. What this means is that (a) and (b) become the same and anything is included in the background set. This will be a very big background. 
sample_expr_lim_cancer=2 #0 
# Number passed to --sample-expr-support-cancer 
# Above we have the number of samples required in the cancer target sample. This means that a kmer is included in the foreground set if its (normalized) expression is >=2 reads. 
cohort_expr_lim_cancer=0  #'1' '5'
# Number passed to --cohort-expr-support-cancer
# Above is the expression value required in at least n samples of the cancer cohort. This is a recurrence parameter
expr_n_limit_cancer=1 #'2' '10' 'none'
# Number passed to --n-samples-lim-cancer
# Above we are specifying the n samples in the cancer cohort. 
#For example --cohort-expr-support-cancer=0 and --expr_n_limit_cancer 2 means that a kmer is included in the foreground set if it is expressed with >0 reads (strict inequality because zero) in >=2 cancer cohort samples 


### Management tip: Intermediate running stage 
# What will happen is as follow: (Times below are given for a background of GTEXcore on the full set of genes) 
# 1) The background is loaded and filtered with --cohort-expr-support-norm --n-samples-lim-normal. The output is written to disk as an intermediate file (1). The intermediate file (1) is pretty long to compute in the real case. 5 days or so. 
# 2) The foreground is loaded and filtered according to --cohort-expr-support-cancer --expr_n_limit_cancer. The output is written to disk as an intermediate file (2). The intermediate file (2) is not too long to compute 1h to 1 day. 
# 3) The rest of the computation, will filter the foreground according to --sample_expr_lim_cancer, output some summary files (count files), and perform the operation foregroud minus background. The kmers remaining after this step are saved to disk as final result. This is quickly computed 4 hours. 
# The intermediate stages allow you to speed up the computation. 
# It is suboptimal to runs many combinations of --cohort-expr-support-norm --n-samples-lim-normal. --cohort-expr-support-cancer --expr_n_limit_cancer. --sample_expr_lim_cancer from the start. Because it will just recompute the same Intermediate file (1) several times. And this is slow to compute. 
# The optimal way is to run: 
# a) one paramater pair for --cohort-expr-support-norm --n-samples-lim-normal. Fix any value of interest for the rest. Wait until it finishes. The intermediate file (1) will be computed. One intermediate file (2) will be computed. And one result file will be produced for the given parameters setting. 
# b)  Then take the same parameter pair for --cohort-expr-support-norm --n-samples-lim-normal. And start playing with many values for --cohort-expr-support-cancer --expr_n_limit_cancer. Fix any value for the rest (--sample_expr_lim_cancer). Wait until it finishes. This will reuse the intermediate file (1), so we skip the long computation for (1). This will create the intermediate files (2) for all the values of interest. And will create the output files for the given parameters combinations. 
# c) Finally we take the same parameters values for --cohort-expr-support-norm --n-samples-lim-normal, and the parameters values for --cohort-expr-support-cancer --expr_n_limit_cancer used before. Start playing with the value of --expr_n_limit_cancer. This will reuse the intermediate file (1), so we skip the long computation for (1). This will reuse the intermediate files (2), so we skip computations again. This will only compute the step 3) and will produce an output for the given parameter settings. 

### Runs 
while read sample; do
	for mutation_canc in ref; do #We will limit ourselves to the cancer case in this project 
		## Organize folders
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
		logfile=${log_dir}/${sample}.cancerspec.${suffix}.${tag_normals}.log
		## Extract parameters	
		cohort_expr_lim_normal=$(echo $expr_nsamples_limit_normal | cut -f1 -d ',')
		expr_n_limit_normal=$(echo $expr_nsamples_limit_normal | cut -f2 -d ',')
		
		if [ ${expr_n_limit_cancer} == 'none' ]; then
			cohort_expr_lim_cancer='none'
		fi
		test_output_exist=$(echo ${output_dir}/G_${sample}_${mutation_canc}_SampleLim${sample_expr_lim_cancer}.0CohortLim${cohort_expr_lim_cancer}.0Across${expr_n_limit_cancer}_FiltNormals${tag_normals}Cohortlim${cohort_expr_lim_normal}.0Across${expr_n_limit_normal}_FiltUniprot.tsv | sed 's,Any,None,g' | sed 's,none,None,g' |sed 's/None\.0/None/g' )
		
		## Cmd
		cmd000="immunopepper cancerspecif --cores $parallel --mem-per-core $mem --kmer $kmer --expression-fields-c "segmentExpr" "junctionExpr" --path-cancer-matrix-edge ${input_Junc_cancer} --ids-cancer-samples "${sample}" --mut-cancer-samples ${mutation_canc} --whitelist-cancer ${whitelist_cancer} --path-cancer-libsize ${libsize_cancer} --normalizer-cancer-libsize ${normalizer_cancer_libsize} --whitelist-normal ${whitelist_normal} --path-normal-libsize ${libsize_normal} --normalizer-normal-libsize ${normalizer_normal_libsize} --output-dir $output_dir --sample-expr-support-cancer ${sample_expr_lim_cancer} --uniprot ${uniprot} --parallelism ${parallelism} --out-partitions ${out_partitions} --path-normal-matrix-segm ${input_Segm_normal} --path-normal-matrix-edge ${input_Junc_normal} --path-normal-kmer-list ${input_annot_cancer} --output-count ${file_count} --interm-dir-norm ${output_norm} --interm-dir-canc ${output_canc} --tag-prefix 'G'" #TODO add back scratch for cancer? --scratch-dir 'TMPDIR'" #TODO output count remove? 
		
		## Here we set parameters that are not systematically included depending on the user's request.  
		if [ ${cohort_expr_lim_normal} != 'Any' ]; then
			cmd00="${cmd000} --cohort-expr-support-norm ${cohort_expr_lim_normal}"
		else
			cmd00="${cmd000}"
		fi
		if [ ${expr_n_limit_normal} != 'Any' ]; then
			cmd0="${cmd00} --n-samples-lim-normal ${expr_n_limit_normal}"
		else
			cmd0="${cmd00}"
		fi
		
		if [ ${expr_n_limit_cancer} != 'none' ]; then 	
			cmd1="${cmd0} --cohort-expr-support-cancer ${cohort_expr_lim_cancer} --n-samples-lim-cancer ${expr_n_limit_cancer}"
		else 
			cmd1="${cmd0}"
		fi 
		
		## Here we redirect the printing of the output 
		cmd3="${cmd1} > ${output_dir}/${sample}.${mutation_canc}.run_cancerspecif.${suffix}.Cec${cohort_expr_lim_cancer}.Cn${expr_n_limit_cancer}.Ces.${sample_expr_lim_cancer}.Nec${cohort_expr_lim_normal}.Nn${expr_n_limit_normal}_nbtc.log 2>&1"	
		
		## Send runs matching conditions. I skipped anything related to the "batch" feature as it is complex and not needed here. Batch runnig allows to send the runs in several batches for big files	
		## Run
		if [ "$local_" = "print_only_the_command" ] ; then
			echo $cmd3
		else
			if [ ! -f "${test_output_exist}/_SUCCESS" ] ; then
				echo $test_output_exist	
				echo $cmd3
				echo $cmd3 | bsubio -n ${parallel} -J ${sample_back} -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -R "rusage[scratch=$scratch_mem]" -o $logfile #-e ${logfile}.e -o $logfile #-R "span[hosts=1]" -o $logfile
			fi
		fi
	done
done < ./tmp_samples #/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/BRCA_5samples_spladder_full.csv #./tmp_samples
