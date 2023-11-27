#!/bin/bash

set -e

mem=6000
basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results

junc_dirs="${basedir}/junctions_spladder_filtered ${basedir}/junctions_spladder"
for junc_dir in $junc_dirs
do
    outdir=${junc_dir}_projected
    mkdir -p $outdir
    for chr in $(seq 1 22) X Y M
    do
        outfname=${outdir}/$(basename $junc_dir).projected.chr${chr}.hdf5
        if [ ! -f ${outfname} ]
        then
            echo "python $(pwd)/collect_single_intron_hdf5.py ${junc_dir} chr${chr} ${junc_dir}.all_coords.sorted.uniq.tsv.gz ${outfname}" | bsub -n 1 -J gtex_prj -M ${mem} -o /dev/null -We 8:00 -R "rusage[mem=${mem}]"
        fi
    done
done
