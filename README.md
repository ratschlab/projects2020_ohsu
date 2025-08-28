This is the analysis code of the preprint from Prelot and David et. al.

## PEPTIDE GENERATION STEPS
Generation of foreground peptides: 
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_generation/foreground
Generation of GTEX peptides
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_generation/background

## FILTERING STEPS 
Codes in https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/python_pipeline 
Steps explained in https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/python_pipeline/README.txt
Note: This is heavily inspired from https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/spark_pipeline, but this pipeline has not been used to generate the final results

## FASTA GENERATION #*pool*.fa
Codes in
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/20230523_speedup_kmer_to_metadata_matching.ipynb
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/p20230320_format_kmers_to_peptides_FASTA-may2023.ipynb

## EXPERIMENT PER PEPTIDES #*experiments_per_peptide*
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/p20230503_parse_fasta_to_experiments_per_peptide.ipynb
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/helpers_map.py
(The old development version is in
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/debug/single_case/p20230726_extract_exp_map_for_fasta.ipynb)

## EXPERIMENT MAPS *experiment_map*
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptide_search_format/notebooks/p20230726_parse_experiment_map.ipynb

## COUNT FILES #*_samp_chrt_norm_mot_unip.tsv*
Generated along the
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/peptides_filtering/python_pipeline
Added the star filtering to the files here
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF/20240204_count_files-for-julianne-figures.ipynb

## PROTEOMICS STEPS 
Codes in https://github.com/ratschlab/projects2020_ohsu/tree/master/pnnl
Steps explained in https://github.com/ratschlab/projects2020_ohsu/tree/master/pnnl/README_proteomics.txt

## PLOTTING

### Parse the filtering results
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/p20230315_extract_filtering_results.ipynb
### Plot the filtering results (intersect the filtered OHSU and ETH kmers on a barplot)
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/20240320_barplot_Experiment_intersection-filtered-kmers-plottingOnly.ipynb

### Parse the proteomics results
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/send_e-p-k.sh
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/run_extract-peptides-to-kmer.sh
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/extract-peptides-to-kmer.py
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/run_example.sh

### Plot the proteomics results (Intersect validated OHSU and ETH kmers on a barplot)
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/20240116_barplot_Experiment_intersection-peptides-to-kmer-plottingOnly.ipynb

### Helpers in 
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/plotting/OHSU_collab_Rerun_allRF_proteomics/

### Perform the intersection of the generated peptides for each pipeline 
https://github.com/ratschlab/projects2020_ohsu/tree/master/eth/MY_Master_thesis_rerun_LP

