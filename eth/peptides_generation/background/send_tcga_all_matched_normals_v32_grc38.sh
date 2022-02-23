#!bash/bin 

mode_send='annot' #'all'
job_name=TCGA_all_normals #TCGA_BRCA_normals  #TCGA_all_normals
cap=1000 
target=mx_commit_438849e_pya.0.17.1_ref_mode_cap${cap}
local_=non_local  #$3 # "run_local"
parallel=4
batch_size=10 
time_=120
mem=8000
batch=10

basedir=/cluster/work/grlab/projects/projects2020_OHSU
log_dir=./logs_${target}
mkdir -p ${log_dir}

### Immunopepper Run
annotation=${basedir}/annotation/gencode.v32.annotation.gtf
genome=${basedir}/genome/GRCh38.p13.genome.fa
if [ ${job_name} == 'TCGA_all_normals' ]; then 
    splice_path=${basedir}/spladder/TCGA_all_normals/genes_graph_conf2.merge_graphs.pickle
    count_path=${basedir}/spladder/TCGA_all_normals/genes_graph_conf2.merge_graphs.count.hdf5
elif [ ${job_name} == 'TCGA_BRCA_normals' ]; then 
    splice_path=/cluster/work/grlab/projects/projects2020_OHSU/spladder/TCGA_BRCA_normals/spladder_confidence_2/genes_graph_conf2.merge_graphs.pickle
    count_path=/cluster/work/grlab/projects/projects2020_OHSU/spladder/TCGA_BRCA_normals/spladder_confidence_2/genes_graph_conf2.count.hdf5
fi
outdir=/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/${job_name}/${target}
mkdir -p ${outdir}
mutation=ref

for sample in TCGA-06-0678-11A.all ; do
base_cmd="immunopepper build --samples ${sample} --output-dir $outdir --ann-path ${annotation} --splice-path ${splice_path} --count-path ${count_path} --ref-path ${genome} --kmer 9 --mutation-mode ${mutation} --verbose 1 --parallel ${parallel} --batch-size ${batch} --complexity-cap $cap "

if [ ${mode_send} == 'all' ]; then  
    cmd1="${base_cmd} --process-chr "chr22"  --use-mut-pickle --all-read-frames --cross-graph-exp " #Note: no libsize custom for this case 
elif [ ${mode_send} == 'annot' ]; then 
    cmd1="${base_cmd}  --use-mut-pickle  --cross-graph-exp"
fi 

cmd="${cmd1} > ${outdir}/${sample}_run_peptides.${mutation}.log 2>&1"

echo $cmd
#echo $cmd | bsub -G ms_raets -n ${parallel} -J ${job_name} -W ${time_}:00 -R "rusage[mem=${mem}]" -R "span[hosts=1]" -o ${log_dir}/${job_name}_${sample}_run_peptides.${mutation}.lsf.log

done
