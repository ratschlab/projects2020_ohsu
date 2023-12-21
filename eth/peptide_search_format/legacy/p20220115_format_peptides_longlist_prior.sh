sample=TCGA-AO-A0JM-01A-21R-A056-07

base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102
cd $base_dir

# Generate all the longlists ${sample}.all_d4aee54_GTEXcore_kmer_longlist.tsv
echo 'longlist'
for folder in filter*; do echo $folder; cat $folder/commit_d4aee54_GTEXcore/*Uniprot*/*part* | cut -d '	'  -f1 | sort | uniq  >  $folder/commit_d4aee54_GTEXcore/$( echo ${folder} | sed 's,filter_,,g')_d4aee54_GTEXcore_kmer_longlist.tsv ; done

# Generate the metadata for the sample of interest
for sample in TCGA-BH-A18V-01A-11R-A12D-07 TCGA-C8-A12P-01A-11R-A115-07 TCGA-A2-A0D2-01A-21R-A034-07 TCGA-A2-A0SX-01A-12R-A084-07 TCGA-AO-A0JM-01A-21R-A056-07; do 
echo "peptide metadata for sample ${sample}"
for kmer in $(cut -f1 -d '	' filter_${sample}.all/commit_d4aee54_GTEXcore/${sample}.all_d4aee54_GTEXcore_kmer_longlist.tsv); do echo ${kmer},$(grep $kmer cohort_mutNone/meta_peptide_pooled_pq.tsv) >> filter_${sample}.all/commit_d4aee54_GTEXcore/${sample}.all_d4aee54_GTEXcore_kmer_peptides_raw.tsv; done
done 


# Temporary fix Generate the metadata for the sample of interest
for sample in TCGA-BH-A18V-01A-11R-A12D-07 TCGA-C8-A12P-01A-11R-A115-07 TCGA-A2-A0D2-01A-21R-A034-07 TCGA-A2-A0SX-01A-12R-A084-07 TCGA-AO-A0JM-01A-21R-A056-07; do 
echo "peptide metadata for sample ${sample}"
for kmer in $(cut -f1 -d '	' filter_${sample}.all/commit_d4aee54_GTEXcore/${sample}.all_d4aee54_GTEXcore_kmer_longlist.tsv); do echo ${kmer},$(grep $kmer /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v3_d2d2574_conf2_annotFrame_cap0_runs/TCGA_Breast_1102/cohort_mutNone/meta_peptide_pooled_pq.tsv) >> filter_${sample}.all/commit_d4aee54_GTEXcore/${sample}.all_d4aee54_GTEXcore_kmer_peptides_raw_flag.tsv; done
done 
