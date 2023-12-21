### Documentation of the scripts used for the peptide neighbor search
Note: There are 3 experimental setup
1) Union of all experimental conditions per sample and per pipeline 
- "Pool" ETH
- "Pool" OHSU
2) Union of all experimental conditions per sample and across 2 pipelines
- "Pool" ETH UNION "Pool" OHSU => ETH "Union" OHSU 
3) Single experimental conditions per sample and per pipeline

## Code Detail
1) Need to extract the tryptic peptides from Fasta files
./20230712-extract-peptides:
runall.sh
script_suite.sh

2) Index the irrelevant, index the relevants, search the neighbors of relevant
in irrelevant, concatenate relevants and neighbors
./20230713-tide-index:
runall-createIndex-eth-ohsu.sh
script_index.sh

3) Search the relevant and neighbors
./20230714-tide-search:
script_search.sh
tide-search-rewriten.sh

4) Extract the PSM from single experiments
./20231212-FDR-correct:
psm_extract

4.bis) Extract the PSM from "Pool" and "Union" experiments
./20231212-FDR-correct:
vanilla_extract

5) Run the FDR and assign confidence 
./20231212-FDR-correct
conf_assign

6) Extract the PSM from each experiments from "Pool" and "Union" files
TODO

## Other relevant codes
andy_bin:
extract-peptides.py #Trypsine
filter-db-search.py # Search
filter-peptides.py # Trypsine
filter-similar-peptides.py # Index
pep-dist.py
plot-pep-cdf.py
rerank-filter-search.py # Unused, reimplemented in
./20231212-FDR-correct/psm_extract 
subset-neighbor-search # Neighbors

