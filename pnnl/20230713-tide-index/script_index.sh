#!/bin/bash

#SBATCH --job-name=index
#SBATCH --output=index.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2GB

crux_home=$1
scripts_home=$2
folder_eth=$3
folder_ohsu=$4
overwrite=$5
h_sapiens_fasta=$6

ppmError=40

# generate irrelevant peptides
${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/irrelevant ${h_sapiens_fasta} tide-indicies/irrelevant
echo

# generate relevant index
# TODO NOTE USING -filter-unique.fasta here and -filter.fasta below.
# This is done to generate neighbors faster
cat ${folder_eth}/peptide-extracted-filter-unique.fasta ${folder_ohsu}/peptide-extracted-filter-unique.fasta >relevant-peptides.fasta
${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/relevant relevant-peptides.fasta tide-indicies/relevant 
echo

# generate neighbors
python3 ${scripts_home}/subset-neighbor-search/pepsim.py --min-score 0.25 --frag-bin-size 0.02 --mz-thresh $ppmError tide-indicies/relevant/tide-index.peptides.txt tide-indicies/irrelevant/tide-index.peptides.txt >simPeptides.txt
echo

# legacy command code. In practice this does not do anything except convert to a fasta file
python3 ${scripts_home}/filter-similar-peptides.py simPeptides.txt 0.25

# combine neighbors and relevant
# TODO same as note above
cat simPeptides_filter.fa relevant-peptides.fasta >finalDb.fasta
echo

${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/final finalDb.fasta tide-indicies/final

