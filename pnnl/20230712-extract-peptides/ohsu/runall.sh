#!/bin/bash

# folder containing ETH kmer fasta 
fol=/eternity/data/eth_pnnl_ohsu_datashare/andy/results/20230712-eth-ohsu-brca/data/ohsu

for f in $fol/*fa;
do
	echo $f
	curF=$(basename $f _pool_kmer.fa)
	mkdir $curF
	cd $curF

	python3 ~/projects/ohsu-eth-collab/bin/extract-peptides.py $f >log.txt
	python3 ~/projects/ohsu-eth-collab/bin/filter-peptides.py peptide-extracted.fasta >>log.txt
	python3 ~/projects/ohsu-eth-collab/bin/plot-pep-cdf.py peptide-extracted-filter.fasta >>log.txt

	cd ../
done
