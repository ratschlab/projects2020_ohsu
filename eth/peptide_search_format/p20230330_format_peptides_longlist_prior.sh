
basedir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/

# Generate all the longlists #all samples?
cd ${basedir}filtering_samples/filters_22March_order_wany_wAnnot
zcat G_TCGA*Filt*tsv.gz |  cut -f1-2 -d '       ' |grep -v kmer | smt_peptides_longlist_prior.shrt | uniq > G_TCGA_Allsamples_all_exp_kmer_longlist.tsv

# Generate the metadata for the samples of interest
cd ${basedir}
zcat cohort_mutNone/tmp_out_ref_batch*/ref_sample_peptides_meta.gz > tmp_pooled_ref_sample_peptides_meta.gz
rm G_TCGA_Allsamples_all_exp_kmer_metadata_raw.tsv; while read line; do kmer=$(echo ${line} | cut -f1 -d ' ' ); coord=$(echo ${line}| cut -f2 -d ' ') ; echo "${kmer}\t${coord}\t$(grep $kmer ${basedir}/tmp_pooled_ref_sample_peptides_meta.tsv)" >> G_TCGA_Allsamples_all_exp_kmer_metadata_raw.tsv; done < G_TCGA_Allsamples_all_exp_kmer_longlist.tsv
