#!/bin/bash

set -e

mem=15000

#basedir=/cluster/work/grlab/projects/GTEx/rna/results
#basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results
#basedir=/cluster/work/grlab/projects/ICGC/alignments
basedir=/cluster/work/grlab/projects/ICGC/alignments/alignments_gtex_2015-03-29_sparse

### unfiltered case
allcoords=${basedir}.junctions.all_coords.sorted.uniq.tsv.gz
for fname in ${basedir}/*.tsv.gz
do
    outname=${fname%.tsv.gz}.projected.hdf5
    logname=${outname%.hdf5}.log
    if [ ! -f ${outname} ]
    then
        #echo "python $(pwd)/project_intron_tsv_to_single_hdf5.py $fname $allcoords" | bsub -n 1 -J gtex_prj -M $mem -o $logname -We 4:00 -R "rusage[mem=${mem}]"
        sbatch --cpus-per-task 1 --job-name gtex_prj --mem=$mem --output=$logname --time=4:00:00 --wrap "python $(pwd)/project_intron_tsv_to_single_hdf5.py $fname $allcoords"
        #echo "python $(pwd)/project_intron_tsv_to_single_hdf5.py $fname $allcoords"
    else
        echo $fname complete
    fi
done
