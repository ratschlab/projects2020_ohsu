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
Note: For the "UNION" of OHSU and ETH pipelines, a reindexing of the peptides
is done to avoid duplication of the intersection
./20230713-tide-index:
runall-createIndex-eth-ohsu.sh
script_index.sh

3) Search the relevant and neighbors
./20230714-tide-search:
script_search.sh
tide-search-rewriten.sh

4) Extract the PSM and compute the FDR 
./20231212-FDR-correct:
4.1) Extract the PSM for the "POOL" and the "Union"
 ./20231212-FDR-correct/POOL_1_concat_extract/
 run_all_Vextract.sh
 script_extract.sh
4.2) Run the FDR for the "POOL" and the "Union" 
 ./20231212-FDR-correct/POOL_2_conf_assign
 runall_conf_union_pool.sh
 script_confidence.sh
4.3) Only for the "UNION" case: Separate peptides from OHSU and ETH. Reindex
./20231212-FDR-correct/POOL_3_FDR_divide
20240220_debug_joint_proteomics_split.py
run_example.sh
runall_divide.sh
script_divide.sh
4.4) Separate the hits per experiment
./20231212-FDR-correct/POOL_3_FDR_extract
FDR_to_experiments.py
helpers_split.py
runall_split.sh
script_split.sh
4.a) Extract the PSM for the "single" experiments (Need to re-rank and
separate per experiment)
./20231212-FDR-correct/SINGLE_1_psm_extract
helpers_psm.py
psm_to_experiments.py
runall_psm.sh
script_psm.sh
4.b) Run the FDR for every single experiment
./20231212-FDR-correct/SINGLE_2_conf_assign
runall_conf_single.sh
script_confidence.sh


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

