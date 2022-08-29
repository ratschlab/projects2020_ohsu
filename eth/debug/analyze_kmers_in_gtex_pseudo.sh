# Extract metadata for a list of problematic kmers (brca)
# cd /cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/TCGA_Breast_1102/cohort_mutNone
#for line in $(cut -f2 '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/TCGA_Breast_1102/filter_TCGA-BH-A18V-01A-11R-A12D-07.all/commit_372a147_full_flags_GTEXcore/recurrence>100_expression>10.tsv') ; do grep $line ref_sample_peptides_meta.csv >../test_debug/${line}_peptides.txt ;  done

# Get the batch to gene equivalence (gtex) 
# In python 
# >>> data, meta = pickle.load(open('/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results/splicing/spladder/genes_graph_conf2.merge_graphs.pickle', 'rb'))
#>>> gene_to_batch = [[i, j.name]   for i, j in enumerate(data) ]
#>>> bar = pd.DataFrame(gene_to_batch)
#bar.to_csv('/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/batch_to_gene.txt',  index = None)



# File extracts from the metadatafile. One per kmer. 
run_path=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/commit_v3_TEST_merged3_372a147_medium_run_conf2_annotFrame_cap0_runs/TCGA_Breast_1102/test_debug

cd $run_path

batchtogene=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/batch_to_gene.txt

gtex_base=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/cohort_mutNone/tmp_out_ref_batch_



# Option 1 get batches
for file_ in *peptides.txt; do if [ -s ${file_} ] ; then  gene=$(head -1 ${file_} | cut -f5 -d ',' ) ; batch=$(grep ${gene} ${batchtogene} | cut -f1 -d ',') ; echo $batch; fi; done

#Option3
for file_ in *peptides.txt; do if [ -s ${file_} ] ; then kmer=$(echo ${file_} | cut -f1 -d '_');  gene=$(head -1 ${file_} | cut -f5 -d ',' ) ; batch=$(grep ${gene} ${batchtogene} | cut -f1 -d ',') ; cd ${gtex_base}$batch; parquet-tools csv ref_graph_kmer_JuncExpr.pq | grep ${kmer} > tmp_${kmer}_JuncExpr.csv; parquet-tools csv ref_graph_kmer_SegmExpr.pq | grep ${kmer} > tmp_${kmer}_SegmExpr.csv; parquet-tools csv ref_sample_peptides_meta.pq | grep ${kmer} > tmp_${kmer}_sample_peptides_meta.pq; cd $run_path; fi; done
