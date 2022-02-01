

base_dir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/v2_v2.5f0752a_conf2_annotFrame_cap0_runs_pya0.17.1/TCGA_Breast_1102
cd $base_dir

for folder in filter*; do echo $folder; cat $folder/commit_d4aee54_GTEXcore/*Uniprot*/*part* | cut -d '     ' -f1 | sort | uniq  >  $folder/commit_d4aee54_GTEXcore/$( echo ${folder} | sed 's,filter_,,g')_d4aee54_GTEXcore_kmer_longlist.tsv ; done

for kmer in $(cut -f1 -d '     ' filter_TCGA-A2-A0D2-01A-21R-A034-07.all/commit_d4aee54_GTEXcore/TCGA-A2-A0D2-01A-21R-A034-07.all_d4aee54_GTEXcore_kmer_longlist.tsv); do echo ${kmer},$(grep $kmer cohort_mutNone/meta_peptide_pooled_pq.tsv) >> filter_TCGA-A2-A0D2-01A-21R-A034-07.all/commit_d4aee54_GTEXcore/TCGA-A2-A0D2-01A-21R-A034-07.all_d4aee54_GTEXcore_kmer_peptides_raw.tsv; done
