#!/bin/bash
sample_type="brca"
STAR='_STAR'
filter_pattern=filtering_samples/filters_19May_order_5ge_wAnnot
make_longlist=False
extract_metadata=True
if [ "${sample_type}" == "ov" ] ; then 
echo 'ov'
base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374
sample_list="TCGA-25-1319-01A-01R-1565-13 TCGA-25-1313-01A-01R-1565-13 TCGA-61-2008-01A-02R-1568-13 TCGA-24-1431-01A-01R-1566-13 TCGA-24-2298-01A-01R-1569-13"
elif [ "${sample_type}" == "brca" ] ; then
echo 'brca'
base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102
sample_list="TCGA-AO-A0JM-01A-21R-A056-07 TCGA-BH-A18V-01A-11R-A12D-07 TCGA-C8-A12P-01A-11R-A115-07 TCGA-A2-A0SX-01A-12R-A084-07 TCGA-A2-A0D2-01A-21R-A034-07"
fi

### Main 

if [ "${make_longlist}" == "True" ]; then 
	cd $base_dir
	for sample in ${sample_list}; do  
		# Generate all the longlists ${sample}_kmer_longlist.tsv.gz
		echo 'longlist'
		for folder in ${filter_pattern}*; do echo $folder; zcat $folder/*${sample}*FiltNormals*tsv.gz | cut -d '	'  -f1-2 |grep -v kmer | sort | uniq  >  $folder/G_${sample}_kmer_longlist.tsv ; gzip  $folder/G_${sample}_kmer_longlist.tsv ; done
	done
fi 

if [ "${extract_metadata}" == "True" ]; then
	for sample in ${sample_list}; do
		metadata_pool=${base_dir}/tmp_pooled_ref_sample_peptides_meta.tsv
		extracted_metadata=${base_dir}/${filter_pattern}/G_${sample}_grep_metadata_raw.tsv
		input_longlist=${base_dir}/${filter_pattern}${STAR}/G_${sample}_kmer_longlist.tsv.gz
		input_longlist2=${base_dir}/${filter_pattern}${STAR}/G_${sample}_kmer_longlist.tsv
		rm ${extracted_metadata}
		echo "...extract metadata ${filter_pattern}${STAR} $sample"
		gunzip ${input_longlist}
		while read line; do kmer=$(echo ${line} | cut -f1 -d ' ' ); coord=$(echo ${line}| cut -f2 -d ' ') ; echo "${kmer}\t${coord}\t$(grep $kmer ${metadata_pool})" >> ${extracted_metadata}; done < ${input_longlist2}
		gzip ${extracted_metadata}
		echo "${extracted_metadata}.gz"
		gzip ${input_longlist2}
	done 
fi
