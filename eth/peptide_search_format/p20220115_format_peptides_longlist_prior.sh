sample=TCGA-AO-A0JM-01A-21R-A056-07

base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102
cd $base_dir

# Generate all the longlists ${sample}.all_d4aee54_GTEXcore_kmer_longlist.tsv
echo 'longlist'
#for folder in filter*; do echo $folder; cat $folder/commit_d4aee54_GTEXcore/*Uniprot*/*part* | cut -d '	'  -f1 | sort | uniq  >  $folder/commit_d4aee54_GTEXcore/$( echo ${folder} | sed 's,filter_,,g')_d4aee54_GTEXcore_kmer_longlist.tsv ; done

# Generate the metadata for the sample of interest
echo "peptide metadata for sample ${sample}"
for kmer in $(cut -f1 -d '	' filter_${sample}.all/commit_d4aee54_GTEXcore/${sample}.all_d4aee54_GTEXcore_kmer_longlist.tsv); do echo ${kmer},$(grep $kmer cohort_mutNone/meta_peptide_pooled_pq.tsv) >> filter_${sample}.all/commit_d4aee54_GTEXcore/${sample}.all_d4aee54_GTEXcore_kmer_peptides_raw.tsv; done
