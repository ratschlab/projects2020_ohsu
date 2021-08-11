#!/bin/bash
set -e

mem=10000
time_=24

cap=0 #TODO change 
target=juliannetestbug_pya0.17.1_annot_cap${cap}_p_spladder2
local_=run_local  #$3 # "run_local"
parallel=4 $2
batch_size=10 $4


basedir=/cluster/work/grlab/projects/projects2020_OHSU
base_path=${basedir}/peptides_generation
log_dir=./logs_${target}
mkdir -p ${log_dir}


### Immunopepper Run 
whitelist_genes=/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/genes_coding_gencode_v32_inter_v30_wo_version.txt
annotation=/cluster/work/grlab/projects/projects2020_OHSU/annotation/gencode.v32.annotation.gtf #/cluster/work/grlab/projects/projects2020_OHSU/annotation/gencode.v32_IntersectGenesInV30.gtf
genome=${basedir}/genome/GRCh38.p13.genome.fa
vcf_path="${basedir}/germline_variants/mergedfiles_clean_stringentfilter.matchIds.h5" # Dummy, see if we have the variant calls with the right genome 
maf_path="${basedir}/somatic_variants/pancan.merged.v0.2.6.PUBLIC.matchIds.maf" # Dummy see if we have the variant call with the same geno,e
heter_code='0'
kmer='9'


#WARNING: using whitelist genes and libsize custom
#source deactivate
#source activate myimmuno3
echo "WARNING check activation myimmuno3"

for cancer_type in OV ; do # TODO add back the OV #OV; do 
	for cf_level in 1 ; do #TODO add #2; do
	       for read_frame in annot ; do # TODO add all; do 	
			while read sample; do
				outdir=${base_path}/${target}/TCGA_${cancer_type}_cancers/spladder_confidence_${cf_level}/${read_frame}_frames
                       		mkdir -p $outdir
				for mutation in ref; do #somatic_and_germline  germline somatic ref ; do  
				out_1=${outdir}/${mutation}_sample_${kmer}mer.pq.gz
                		out_2=${outdir}/${mutation}_annot_${kmer}mer.pq.gz
                		out_3=${outdir}/${mutation}_annot_peptides.fa.pq.gz
                	        out_4=${outdir}/${mutation}_sample_peptides.fa.pq.gz
                	        out_5=${outdir}/${sample}/${mutation}_sample_peptides_meta.tsv.gz.pq
                		if [ ! -f $out_1 ] || [ ! -f $out_2 ] || [[ ! -f $out_5 ]]; then 
					
					echo "${outdir}/${mutation}${out_file} does not exist"
				       #TODO REPLACE THIS WITH THE MERGED GRAPH 
				       # TODO add the softlink path 	spladder_path=${basedir}/spladder/TCGA_${cancer_type}_cancers/spladder_confidence_${cf_level}/genes_graph_conf${cf_level}.${sample}
				       spladder_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Ovarian_374_conf${cf_level}_results/splicing/spladder/genes_graph_conf${cf_level}.${sample}.pickle
				       count_path=/cluster/work/grlab/projects/projects2021-immuno_peptides/results/TCGA_for_neoepitopes/TCGA_Ovarian_374_conf${cf_level}_results/splicing/spladder/genes_graph_conf${cf_level}.merge_graphs.${sample}.count.hdf5
				       cmd_base="immunopepper  build --verbose 1 --output-samples ${sample} --output-dir ${outdir} --ann-path ${annotation} --splice-path ${spladder_path} --count-path ${count_path} --ref-path ${genome} --kmer ${kmer} --mutation-mode ${mutation} --somatic ${maf_path} --germline ${vcf_path} --batch-size ${batch_size} --complexity-cap $cap --genes-interest ./tmp_gene" #TODO add whitelist genes --genes-interest ${whitelist_genes}"
				      
				       if [ "$read_frame" == "all" ]; then 
						cmd="${cmd_base} --parallel ${parallel} --process-chr "chr22" --use-mut-pickle --all-read-frames  > ${outdir}/${sample}_run_peptides.${mutation}.log 2>&1"
					else
						cmd="${cmd_base} --parallel ${parallel} --use-mut-pickle > ${outdir}/${sample}_run_peptides.${mutation}.log 2>&1"
					fi



					if [ "$local_" = "run_local" ] ; then
						echo "running local"
						echo $cmd
					else
						echo $cmd
						echo $cmd | bsub -J for_annot -n ${parallel} -J ${mutation}_ohsu -W ${time_}:00 -R "rusage[mem=${mem}]" -o ${log_dir}/${sample}_run_peptides.${mutation}.lsf
					fi
				fi
				done
			done < ./tmp_send_file # TODO add sample file ${basedir}/sample_lists/TCGA_foreground/${cancer_type}_5samples_id_random_final.csv 
		done
	done
done
