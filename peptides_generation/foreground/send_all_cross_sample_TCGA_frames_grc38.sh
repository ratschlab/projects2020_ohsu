#!/bin/bash
set -e

mem=50000
time_=24
local_=run_local #run_cluster 
parallel=1 $2

### Immunopepper Run
cap=1000 #TODO Confirm??
batch_size=10 $4
frames=annot
conf=conf1
basedir=/cluster/work/grlab/projects/projects2020_OHSU
base_path=${basedir}/peptides_generation

if [ "$frame" == "all" ] ; then
	target=v2_a6e5cad_${conf}_allFrame_cap${cap}_runs_pya0.17.1
else
	target=v2_a6e5cad_${conf}_annotFrame_cap${cap}_runs_pya0.17.1
fi
outdir=${base_path}/${target}
log_dir=${outdir}/lsf
mkdir -p $outdir
mkdir -p ${log_dir}


basedir=/cluster/work/grlab/projects/projects2020_OHSU
base_path=${basedir}/peptides_generation
log_dir=./logs_${target}
mkdir -p ${log_dir}

### Immunopepper Run
annotation=${basedir}/annotation/gencode.v32.annotation.gtf
genome=${basedir}/genome/GRCh38.p13.genome.fa
vcf_path="${basedir}/germline_variants/mergedfiles_clean_stringentfilter.matchIds.h5" # Dummy, see if we have the variant calls with the right genome
maf_path="${basedir}/somatic_variants/pancan.merged.v0.2.6.PUBLIC.matchIds.maf" # Dummy see if we have the variant call with the same geno,e
heter_code='0'
kmer='9'

count_path=/cluster/work/grlab/projects/tmp_tinu/TCGA_for_neoepitopes/TCGA_Ovarian_374_conf1_results/splicing_conf1/spladder/genes_graph_${conf}.merge_graphs.count.hdf5
splice_path=/cluster/work/grlab/projects/tmp_tinu/TCGA_for_neoepitopes/TCGA_Ovarian_374_conf1_results/splicing_conf1/spladder/genes_graph_${conf}.merge_graphs.pickle #TODO Link to real path 
sample_file=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_Ov_378.tsv
coding_genes=/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v32.txt

echo "WARNING check activation myimmuno3"

for mutation in ref; do  
	#out_1=${outdir}/${sample}/${mutation}_sample_${kmer}mer.pq
	#out_2=${outdir}/${sample}/${mutation}_annot_${kmer}mer.pq
	#out_3=${outdir}/${sample}/${mutation}_annot_peptides.fa.pq
	#out_5=${outdir}/${sample}/${mutation}_sample_peptides_meta.pq
	#out_4=${outdir}/${sample}/gene_expression_detail.pq

	#if [ ! -f $out_1 ] || [ ! -f $out_2 ] || [[ ! -f $out_5 ]] || [[ ! -f $out_3 ]] || [[ ! -f $out_4 ]]; then 
	while read sample ; do 
	        ## no mutation sample
                if [ "$mutation" == "ref" ]; then
                      sample='cohort'
                fi
		cmd_base="immunopepper  build --verbose 1 --output-samples $(cat ${sample_file} |tr '\n\' '\t') --mutation-sample ${sample} --output-dir ${outdir} --ann-path ${annotation} --splice-path ${splice_path} --count-path ${count_path} --ref-path ${genome} --kmer ${kmer} --mutation-mode ${mutation} --somatic ${maf_path} --germline ${vcf_path} --batch-size ${batch_size} --complexity-cap $cap --genes-interest ${coding_genes}" #TODO Remove tmp genes 
		## Parallel mode
	       if [ "$parallel" -gt 1 ]; then 
			cmd1="${cmd_base} --parallel ${parallel} --use-mut-pickle --cross-graph-expr" 
		else
			cmd1="${cmd_base} --use-mut-pickle --cross-graph-expr"
		fi
		
		## Frame mode
		if [ "$frame" == "all" ] ; then
        		cmd2="${cmd1} --all-read-frames" 
		else
       			cmd2="${cmd1}" 
		fi

		cmd3="${cmd2}  > ${outdir}/mode_build_run_peptides.${mutation}.log 2>&1"

		## Launch 
		if [ "$local_" = "run_local" ] ; then
			echo "running local"
		        echo $cmd3
		else
			echo $cmd3 | bsub -J for_${frame} -n ${parallel} -J ${mutation}_ip_tcga -W ${time_}:00 -R "rusage[mem=${mem}]" -o ${log_dir}/${sample}_run_peptides.${mutation}.lsf 
		fi
	 if [ "$mutation" == "ref" ]; then 
		 break 
	 fi 
      done < ${sample_file}
done
