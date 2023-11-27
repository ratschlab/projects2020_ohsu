#!/bin/bash

set -e

mem=6000

#basedir=/cluster/work/grlab/projects/GTEx/rna/results
basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results

### filtered case
allcoords=${basedir}/junctions_spladder_filtered.all_coords.sorted.uniq.tsv.gz
for fname in ${basedir}/junctions_spladder_filtered/*.tsv.gz
do
    outname=${fname%.tsv.gz}.projected.hdf5
    logname=${outname%.hdf5}.log
    if [ ! -f ${outname} ]
    then
        echo "python $(pwd)/project_intron_tsv_to_single_hdf5.py $fname $allcoords" | bsub -n 1 -J gtex_prj -M $mem -o $logname -We 4:00 -R "rusage[mem=${mem}]"
    else
        echo $fname complete
    fi
done

### unfiltered case
allcoords=${basedir}/junctions_spladder.all_coords.sorted.uniq.tsv.gz
for fname in ${basedir}/junctions_spladder/*.tsv.gz
do
    outname=${fname%.tsv.gz}.projected.hdf5
    logname=${outname%.hdf5}.log
    if [ ! -f ${outname} ]
    then
        echo "python $(pwd)/project_intron_tsv_to_single_hdf5.py $fname $allcoords" | bsub -n 1 -J gtex_prj -M $mem -o $logname -We 4:00 -R "rusage[mem=${mem}]"
    else
        echo $fname complete
    fi
done
