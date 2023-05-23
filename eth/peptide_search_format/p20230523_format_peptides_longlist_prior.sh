#!/bin/bash
sample_type='ov'

if [ "${sample_type} == ov" ] ; then 
base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374
sample_list="TCGA-25-1319-01A-01R-1565-13 TCGA-25-1313-01A-01R-1565-13 TCGA-61-2008-01A-02R-1568-13 TCGA-24-1431-01A-01R-1566-13 TCGA-24-2298-01A-01R-1569-13"
else
base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102
sample_list="TCGA-AO-A0JM-01A-21R-A056-07 TCGA-BH-A18V-01A-11R-A12D-07 TCGA-C8-A12P-01A-11R-A115-07 TCGA-A2-A0SX-01A-12R-A084-07 TCGA-A2-A0D2-01A-21R-A034-07"
fi

cd $base_dir
for sample in ${sample_list}; do  
# Generate all the longlists ${sample}_kmer_longlist.tsv.gz
echo 'longlist'
filter_pattern=filtering_samples/filters_19May_order_5ge_wAnnot
 for folder in ${filter_pattern}*; do echo $folder; zcat $folder/*${sample}*FiltNormals*tsv.gz | cut -d '	'  -f1-2 |grep -v kmer | sort | uniq  >  $folder/G_${sample}_kmer_longlist.tsv ; gzip  $folder/G_${sample}_kmer_longlist.tsv ; done
 
done

