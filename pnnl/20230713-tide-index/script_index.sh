#!/bin/bash

#SBATCH --job-name=index
#SBATCH --output=index.out
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --exclude=compute-biomed-10

crux_home=$1
scripts_home=$2
overwrite=$3
h_sapiens_fasta=$4
dir_irrelevant=$5
folder_pipeline1=$6
folder_pipeline2=$7
union_pipelines=$8

ppmError=40

# generate irrelevant peptides
#${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir ${dir_irrelevant}/tide-indicies/irrelevant ${h_sapiens_fasta} tide-indicies/irrelevant
echo SKIPPING IRRELEVANT INDEX

# generate relevant index
# TODO NOTE USING -filter-unique.fasta here and -filter.fasta below.
# This is done to generate neighbors faster
if [[ ${union_pipelines} == 'T' ]]; then
	python3 ${scripts_home}/20240216_debug_joint_proteomics.py --file-eth ${folder_pipeline1}/peptide-extracted-filter-unique.fasta --file-ohsu ${folder_pipeline2}/peptide-extracted-filter-unique.fasta --save-folder '.' --map-eth-folder ${folder_pipeline1} --map-ohsu-folder ${folder_pipeline2} # Concatenates and re-indexes. Outputs joint-peptide-extracted-filter-unique.fasta
	echo "union pipelines"
	${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/relevant joint-peptide-extracted-filter-unique.fasta tide-indicies/relevant
else
	echo "single pipeline"
	${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/relevant ${folder_pipeline1}/peptide-extracted-filter-unique.fasta tide-indicies/relevant
fi
echo

# generate neighbors
python3 ${scripts_home}/subset-neighbor-search/pepsim.py --min-score 0.25 --frag-bin-size 0.02 --mz-thresh $ppmError tide-indicies/relevant/tide-index.peptides.txt ${dir_irrelevant}/tide-indicies/irrelevant/tide-index.peptides.txt >simPeptides.txt
echo

# legacy command code. In practice this does not do anything except convert to a fasta file
python3 ${scripts_home}/filter-similar-peptides.py simPeptides.txt 0.25

# combine neighbors and relevant
# TODO same as note above
if [[ ${union_pipelines} == 'T' ]]; then
	cat simPeptides_filter.fa relevant-peptides.fasta >finalDb.fasta
else 
	cat simPeptides_filter.fa ${folder_pipeline1}/peptide-extracted-filter-unique.fasta >finalDb.fasta
fi
echo

${crux_home} tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/final finalDb.fasta tide-indicies/final

