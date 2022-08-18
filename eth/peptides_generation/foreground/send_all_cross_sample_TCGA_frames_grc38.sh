#!/bin/bash
set -e

mem=20000
time_=120
local_=run_cluster 
parallel=8 $2

### Immunopepper parameters
start_id=0
cap=0 #TODO 
batch_size=1 $4
frames=annot
conf=conf2
basedir=/cluster/work/grlab/projects/projects2020_OHSU
base_path=${basedir}/peptides_generation
coding_genes=/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/OHSU_gencodev32_proteincodinggeneids.txt
#coding_genes=/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/tmp_genes #TODO update
#coding_genes=./test_genes_recurrent
### Inputs
annotation="${basedir}/annotation/gencode.v32.annotation.gtf"
genome="${basedir}/genome/GRCh38.p13.genome.fa"
vcf_path="${basedir}/germline_variants/mergedfiles_clean_stringentfilter.matchIds.h5" # Dummy, see if we have the variant calls with the right genome
maf_path="${basedir}/somatic_variants/pancan.merged.v0.2.6.PUBLIC.matchIds.maf" # Dummy see if we have the variant call with the same genome
uniprot_kmers="/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv"
heter_code='0'
kmer='9'

#Note 
#/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Breast_1102_results/splicing/spladder/genes_graph_${conf}.merge_graphs.count.hdf5 seems to be the realignement + reurun of spladder graph.
for sample_type in TCGA_Breast_1102; do # TCGA_All_Normals TCGA_Ovarian_374 TCGA_Breast_1102; do 
if [ "$sample_type" == "TCGA_Ovarian_374" ]; then  
    count_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Ovarian_374_results/splicing/spladder/genes_graph_${conf}.merge_graphs.count.rechunked.hdf5
    splice_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Ovarian_374_results/splicing/spladder/genes_graph_${conf}.merge_graphs.pickle  
    sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_Ov_378.tsv
elif [ "$sample_type" == "TCGA_Breast_1102" ]; then
   count_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Breast_1102_results/splicing/spladder/genes_graph_${conf}.merge_graphs.count.hdf5 #RECHUNKED NOT AVAILABLE 
   splice_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Breast_1102_results/splicing/spladder/genes_graph_${conf}.merge_graphs.pickle 
   sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102.tsv
elif [ "$sample_type" == "TCGA_All_Normals" ]; then
   count_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_All_Normals_results/splicing/spladder/genes_graph_${conf}.merge_graphs.count.rechunked.hdf5
   splice_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_All_Normals_results/splicing/spladder/genes_graph_${conf}.merge_graphs.pickle 
   sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_All_Normals/sample_full_All_Normals.tsv
fi

### Outputs
commit=v3_TEST_merged3_57a6a62_libsize #timing_substract_noannot #_libsize #timing_substract
if [ "$frame" == "all" ] ; then
        target=${commit}_${conf}_allFrame_cap${cap}_runs/${sample_type}
else
        target=${commit}_${conf}_annotFrame_cap${cap}_runs/${sample_type}
fi

outdir=${base_path}/${target}
log_dir=${outdir}/lsf
mkdir -p $outdir
mkdir -p ${log_dir}


#TODO remove skip annotation
echo "WARNING check activation myimmuno3"
for mutation in ref; do  
	#out_1=${outdir}/${sample}/${mutation}_sample_${kmer}mer.pq
	#out_2=${outdir}/${sample}/${mutation}_annot_${kmer}mer.pq
	#out_3=${outdir}/${sample}/${mutation}_annot_peptides.fa.pq
	#out_5=${outdir}/${sample}/${mutation}_sample_peptides_meta.pq
	#out_4=${outdir}/${sample}/gene_expression_detail.pq

	#if [ ! -f $out_1 ] || [ ! -f $out_2 ] || [[ ! -f $out_5 ]] || [[ ! -f $out_3 ]] || [[ ! -f $out_4 ]]; then 
	while read sample ; do 
		cmd_base="immunopepper build --verbose 2 --output-dir ${outdir} --ann-path ${annotation} --splice-path ${splice_path} --count-path ${count_path} --ref-path ${genome} --kmer ${kmer}"
	
          
		## Specific processing parameters 
                cmd0="${cmd_base} --cross-graph-expr --skip-tmpfiles-rm --batch-size ${batch_size} --complexity-cap $cap --genes-interest ${coding_genes} --start-id ${start_id} --kmer-database ${uniprot_kmers} --skip-annotation  --libsize-extract" #TODO Remove tmp genes  #Remark, if no output_samples does output all samples from countfile #TODO remove skip annotation

		# mutation mode
		if [ "$mutation" == "ref" ]; then
                      sample='cohort'
		      cmd1="${cmd0}"
		elif [ "$mutation" == "somatic" ]; then
		      cmd1="${cmd0} --mutation-sample ${mutation_sample} --somatic ${maf_path} --use-mut-pickle"
		elif [ "$mutation" == "germline" ]; then  
		      cmd1="${cmd0} --mutation-sample ${mutation_sample} --germline ${vcf_path} --use-mut-pickle"
		fi 

		## Parallel mode
	        if [ "$parallel" -gt 1 ]; then 
			cmd2="${cmd1} --parallel ${parallel}" 
		else
			cmd2="${cmd1}" 
		fi
		

		## Frame mode
		if [ "$frame" == "all" ] ; then
        		cmd3="${cmd2} --all-read-frames" 
		else
       			cmd3="${cmd2}" 
		fi
			
		## Output log 
		cmd_out="${cmd3} > ${outdir}/mode_build_run_peptides.${mutation}.${start_id}.log 2>&1" #--skip-annotation

		## Launch 
		if [ "$local_" = "run_local" ] ; then
			echo "running local"
		        echo $cmd_out
		else
			echo $cmd_out
			echo $cmd_out | bsub -J ohfas${start_id} -n ${parallel} -W ${time_}:00 -R "rusage[mem=${mem}]" -o ${log_dir}/${sample}_run_peptides.${mutation}.${start_id}.lsf 
		fi
	 if [ "$mutation" == "ref" ]; then 
		 break 
	 fi 
      done < ${sample_file}
done
done
