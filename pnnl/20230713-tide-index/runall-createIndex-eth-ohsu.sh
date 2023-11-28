
#!/bin/bash

ppmError=40
overwrite=T
# path for ETH extract peptide folder
folder=/eternity/data/eth_pnnl_ohsu_datashare/andy/results/20230712-eth-ohsu-brca/20230712-extract-peptides

for f in $folder/eth/G*;
do
    echo $f
    curF=$(basename $f)
	tmp_file=${curF/G_/}
    mkdir $tmp_file
    cd $tmp_file

    # generate irrelevant peptides
    ~/projects/ohsu-eth-collab/bin/crux-toolkit/Release/src/crux tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/irrelevant ~/projects/ohsu-eth-collab/fasta/hsapiens.fasta tide-indicies/irrelevant
    echo

    # generate relevant index
    # TODO NOTE USING -filter-unique.fasta here and -filter.fasta below.
    # This is done to generate neighbors faster
	ohsu_file=${curF/G/J}
    cat $folder/eth/$curF/peptide-extracted-filter-unique.fasta $folder/ohsu/$ohsu_file/peptide-extracted-filter-unique.fasta >relevant-peptides.fasta
    ~/projects/ohsu-eth-collab/bin/crux-toolkit/Release/src/crux tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/relevant relevant-peptides.fasta tide-indicies/relevant 
    echo

    # generate neighbors
    python3 ~/projects/ohsu-eth-collab/bin/subset-neighbor-search/pepsim.py --min-score 0.25 --frag-bin-size 0.02 --mz-thresh $ppmError tide-indicies/relevant/tide-index.peptides.txt tide-indicies/irrelevant/tide-index.peptides.txt >simPeptides.txt
    echo

    # legacy command code. In practice this does not do anything except convert to a fasta file
    python3 ~/projects/ohsu-eth-collab/bin/filter-similar-peptides.py simPeptides.txt 0.25

    # combine neighbors and relevant
    # TODO same as note above
    cat simPeptides_filter.fa relevant-peptides.fasta >finalDb.fasta
    echo

    ~/projects/ohsu-eth-collab/bin/crux-toolkit/Release/src/crux tide-index --overwrite $overwrite --mods-spec C+57.02146,K+144.102063 --nterm-peptide-mods-spec X+144.102063 --peptide-list T --output-dir tide-indicies/final finalDb.fasta tide-indicies/final

    cd ../
done
