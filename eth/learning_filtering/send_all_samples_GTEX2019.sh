#!bash/bin 

# This scripts runs immunopepper on the version of GTEX used for the Cancer Cell paper. The genome version used is hg19


### Run parameters
parallel=30
time_=120
mem=20000
batch=1
start_id=30000
job_name=${start_id}_2019
read_frame=annot
mutation_mode=ref #germline


### Immunopepper parameters
#cap=1000
#commit=TEST_empty_loop_282dabc
#commit=v3_TEST_merged3_57a6a62_libsize #substract
#commit=v3_TEST_merged3_372a147_substract
#commit=v3_TEST_merged3_84dc237_substract_withuniprot
#commit=v3_TEST_merged3_372a147_substract_withuniprot
#commit=v3_TEST_merged3_258efee_experim_no_file_exist
#commit=v3_TEST_merged3_258efee_experim_no_file_exist_nouniprot_mem20
#commit=v3_TEST_merged3_372a147_substract_nouniprot_mem50_c10
commit=v3_TEST_merged3_372a147_medium_run
conf=conf2
dir_name=GTEX2019_commit_${commit}_pya.0.17.1_${conf}_${read_frame}_${mutation_mode}_chrall_cap${cap}

gtex_samples=/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_normal_samples_12-2020_nohead.tsv
coding_genes=/cluster/work/grlab/projects/projects2020_OHSU/gene_lists/OHSU_gencodev32_proteincodinggeneids.txt #selection made with JD

outdir=/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/${dir_name}
mkdir -p $outdir
log_dir=${outdir}/lsf
mkdir -p ${log_dir}

heter_code='0'
annotation=/cluster/work/grlab/projects/tumor_profiler/annotation/gencode.v32.annotation.gtf
genome=/cluster/work/grlab/projects/tumor_profiler/genome/GRCh38.p13.genome.fa
uniprot_kmers=/cluster/work/grlab/projects/TCGA/PanCanAtlas/tcga_immuno/uniprot/9mers_uniprot-human-UP000005640_9606.tsv

if [ "${conf}" == "conf2" ] ; then
        # (STAR + SplAdder both use Gencode v32)	
	splice_path=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results/splicing/spladder/genes_graph_conf2.merge_graphs.pickle
	count_path=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results/splicing/spladder/genes_graph_conf2.merge_graphs.count.rechunked.hdf5
fi 

echo $count_path
echo $outdir 

### Runs 
#while read sample; do
out1=${outdir}/${sample}/gene_expression_detail.pq.gz
out2=${outdir}/${sample}/germline_graph_kmer_JuncExpr.pq.gz
out3=${outdir}/${sample}/germline_graph_kmer_SegmExpr.pq.gz
out4=${outdir}/${sample}/germline_sample_peptides_meta.tsv.gz.pq

if [ ! -f $out1 ] || [ ! -f $out2 ] || [ ! -f $out3 ] || [ ! -f $out4 ] ; then 
	cmd="immunopepper build --output-dir $outdir --ann-path ${annotation} --splice-path ${splice_path} --count-path ${count_path} --ref-path ${genome} --kmer "9" --verbose 2 --parallel ${parallel} --batch-size ${batch} --start-id ${start_id} --cross-graph-exp --skip-annotation --skip-tmpfiles-rm --genes-interest ${coding_genes} --kmer-database ${uniprot_kmers}"#--libsize-extract"

    if [ "$read_frame" == 'all' ] ; then 
	    cmd2="${cmd} --all-read-frames"
    else
	    cmd2="${cmd}"
    fi
    cmd3="${cmd2} > ${outdir}/run_peptides.${mutation_mode}_${start_id}.log 2>&1"
   


    echo $cmd3
    #TODO add -R "span[hosts=1]"
    echo $cmd3 | bsub -n ${parallel} -J ${job_name} -W ${time_}:00 -R "rusage[mem=${mem}]" -o ${log_dir}/run_peptides.${mutation_mode}_${start_id}.lsf
fi
#done < $dummy_sample_file
