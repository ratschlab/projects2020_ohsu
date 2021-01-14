#!/bin/bash
set -e

mem=15000
time_=120

target=10c3360_runs_pya0.17.1
local_=non_local  #$3 # "run_local"
parallel=4 $2
batch_size=10 $4


basedir=/cluster/work/grlab/projects/projects2020_OHSU
base_path=${basedir}/peptides_generation
log_dir=./logs_${target}
mkdir -p ${log_dir}


### Immunopepper Run 
annotation=${basedir}/annotation/gencode.v32.annotation.gtf
genome=${basedir}/genome/GRCh38.d1.vd1.phiX.fa
vcf_path="${basedir}/germline_variants/mergedfiles_clean_stringentfilter.matchIds.h5" # Dummy, see if we have the variant calls with the right genome 
maf_path="${basedir}/somatic_variants/pancan.merged.v0.2.6.PUBLIC.matchIds.maf" # Dummy see if we have the variant call with the same geno,e
heter_code='0'
kmer='9'



#source deactivate
#source activate myimmuno3
echo "WARNING check activation myimmuno3"

for cancer_type in BRCA OV; do 
	for cf_level in 1 2; do
	       for read_frame in annot all; do 	
			while read sample; do 
				outdir=${base_path}/${target}/${cancer_type}/${sample}/spladder_confidence_${cf_level}/${read_frame}_frames
                       		mkdir -p $outdir
				for mutation in ref; do #somatic_and_germline  germline somatic ref ; do  
				    if [ ! -f  ${outdir}/${mutation}_junction_${kmer}mer.pq.gz ] || [ ! -f  ${outdir}/${mutation}_back_${kmer}mer.pq.gz ] || [ ! -f  ${outdir}/${mutation}_back_peptides.fa.pq.gz ] ||  [ ! -f  ${outdir}/${mutation}_peptides.fa.pq.gz ] ||  [ ! -f  ${outdir}/${mutation}_metadata.tsv.gz.pq ]; then
				       echo "${outdir}/${mutation}${out_file} does not exist"
				       spladder_path=${basedir}/spladder/${cancer_type}/spladder_confidence_${cf_level}/genes_graph_conf${cf_level}.${sample}
				       cmd_base="immunopepper  build --verbose 1 --samples ${sample} --output-dir ${outdir} --ann-path ${annotation} --splice-path ${spladder_path}.pickle --count-path ${spladder_path}.count.hdf5 --ref-path ${genome} --kmer ${kmer} --mutation-mode ${mutation} --somatic ${maf_path} --germline ${vcf_path} --batch-size ${batch_size}"
				      
				       if [ "$read_frame == "all"" ]; then 
						cmd="${cmd_base} --parallel ${parallel} --use-mut-pickle --all-read-frames  > ${outdir}/${sample}_run_peptides.${mutation}.log 2>&1"
					else
						cmd="${cmd_base} --parallel ${parallel} --use-mut-pickle > ${outdir}/${sample}_run_peptides.${mutation}.log 2>&1"
					fi



					if [ "$local_" = "run_local" ] ; then
						echo "running local"
						$cmd
					else
						echo $cmd
						echo $cmd | bsub -J for_annot -n ${parallel} -J ${mutation}_ohsu -W ${time_}:00 -R "rusage[mem=${mem}]" -o ${log_dir}/${sample}_run_peptides.${mutation}.lsf
					fi
				fi
				done
			done < ${basedir}/sample_lists/TCGA_foreground/${cancer_type}_5samples_id_random_final.csv 
		done
	done
done
