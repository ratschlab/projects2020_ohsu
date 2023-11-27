#!/bin/bash

set -e

basedir=/cluster/work/grlab/projects/GTEx/rna_gencode32_realign/results

### link results - filtered case
echo linking results - filtered cases
mkdir -p ${basedir}/junctions_spladder_filtered
cd ${basedir}/junctions_spladder_filtered
for fname in ${basedir}/alignments/*.conf_2.filt.tsv.gz
do
    ln $fname
done

### link results - unfiltered case
echo linking results - unfiltered cases
mkdir -p ${basedir}/junctions_spladder
cd ${basedir}/junctions_spladder
for fname in ${basedir}/alignments/*.all.tsv.gz
do
    ln $fname
done

### collect all coordinates and make unique
echo collecting coordinates
zcat ${basedir}/alignments/SRR1068687.all.conf_2.filt.tsv.gz | head -n 1 | cut -f -4 | gzip > ${basedir}/junctions_spladder_filtered.all_coords.sorted.uniq.tsv.gz
zcat ${basedir}/alignments/SRR1068687.all.tsv.gz | head -n 1 | cut -f -4 | gzip > ${basedir}/junctions_spladder.all_coords.sorted.uniq.tsv.gz
zcat ${basedir}/alignments/*.conf_2.filt.tsv.gz | grep -v junction_start | cut -f -4 | sort -u | gzip >> ${basedir}/junctions_spladder_filtered.all_coords.sorted.uniq.tsv.gz &
zcat ${basedir}/alignments/*.all.tsv.gz | grep -v junction_start | cut -f -4 | sort -u | gzip >> ${basedir}/junctions_spladder.all_coords.sorted.uniq.tsv.gz &
wait
echo done
