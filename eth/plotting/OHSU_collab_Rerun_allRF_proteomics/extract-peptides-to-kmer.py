import argparse
from Bio import SeqIO
from collections import defaultdict
import glob
import gzip
import numpy as np 
import os 
import pandas as pd
import tarfile

from helpers_barplot_intersection import plot_text, plot_intersection_bars, reader_assign_conf_pep
from helpers_barplot_intersection_kmers import explode_immunopepper_coord, search_result_peptides_ids 
from helpers_barplot_intersection_kmers import get_pep_ids, get_pep_coord, tar_reader
from helpers_barplot_intersection_kmers import validated_filtered_kmers, format_validation_rates
from helpers_barplot_intersection_kmers import compare_OHSU_ETH, kmer_in_bi_exon_peptide


  

def process_proteomics_results(proteomicsdir:str, samples_breast:list, samples_ov:list, \
                               fasta_base_OHSU:str, fasta_base_ETH:str, kmer_files_OHSU:str, \
                               pipelines:list, MS_FDR:str, MS_strategy:str, FDR_limit:float, save_folder:str):

    
    if MS_FDR == '_crema':
        FDR_file = 'crema.peptides.txt'
        col_seq = 'sequence'
        col_qvalue = 'crema q-value'
    elif MS_FDR == '_crux' or  MS_FDR == '':
        FDR_file = 'assign-confidence.target.txt'
        col_seq = 'unmodified sequence'
        col_qvalue = 'tdc q-value'
    else:
        print(f'ERROR: wrong input for {FDR_file}')

    all_samples = []
    all_samples.extend(samples_breast)
    all_samples.extend(samples_ov)

    
    ## Get kmers from result files

    samples_store_kmers = {}
    samples_store_pep = {}
    samples_store_rates_peps = {}
    samples_store_rates_kmers = {}

    for sample in all_samples:
        print('processing', sample)

        sample_short = '-'.join(sample.split('-')[0:3])
        
        samples_store_kmers[sample] = defaultdict(dict)
        samples_store_pep[sample] = defaultdict(dict)
        samples_store_rates_peps[sample] = defaultdict(dict)
        samples_store_rates_kmers[sample] = defaultdict(dict)
        
        for pipeline in pipelines:
            path_single = os.path.join(proteomicsdir, pipeline, sample_short, 
                                       f'assign_conf_per_experiment{MS_FDR}')
            path_pool_pipeline = os.path.join(proteomicsdir, pipeline, sample_short, 
                                              f'assign_conf_pooled_FDR{MS_FDR}')
            path_pool_union = os.path.join(proteomicsdir, 
                                           f'assign_conf_joint_to_{pipeline}{MS_FDR}', sample_short)
            path_TEST_OHSU = os.path.join(proteomicsdir, 'OHSU', sample_short, 
                                       f'assign_conf_per_experiment{MS_FDR}')
            path_TEST_ETH = os.path.join(proteomicsdir, 'ETH', sample_short, 
                                       f'assign_conf_per_experiment{MS_FDR}')

            print(path_single)
            experiment_list = [ i.split('/')[-1] for i in glob.glob(path_single + '/*')] #check
            print(f'Processing {len(experiment_list)} experiments')
            for experiment in experiment_list:
                if pipeline == 'OHSU':
                    original_name = experiment
                    cut_name = experiment[1:]
                else:
                    original_name = experiment
                    cut_name = experiment          

                if os.path.isfile(os.path.join(path_TEST_OHSU, 'J' + cut_name, FDR_file)) and \
                     os.path.isfile(os.path.join(path_TEST_ETH, cut_name, FDR_file)): #Commun experiments               

                    # search 1 experiment, 1 pipeline  
                    if MS_strategy == 'single':
                        df = os.path.join(path_single, original_name, FDR_file)
                    # search all experiments, 1 pipeline
                    if MS_strategy == 'pool':
                        df = os.path.join(path_pool_pipeline, f'tsearch-{original_name}.txt')
                    # search all experiments, 1 union of pipelines
                    if MS_strategy == 'joint':
                        df = os.path.join(path_pool_union, f'tsearch-{original_name}.txt')

                    val, val_rate_tryptic_pep, peptides, df_filtered = reader_assign_conf_pep(df, FDR_limit, col_seq, col_qvalue)

                    df_filtered, val_rate_kmers = validated_filtered_kmers(df_filtered, fasta_base_OHSU, kmer_files_OHSU,
                                                           fasta_base_ETH, sample, experiment, 
                                                           pipeline, col_seq)

                        
                    samples_store_kmers[sample][cut_name][pipeline] = set(df_filtered['kmer'])
                    samples_store_pep[sample][cut_name][pipeline] = peptides
                    samples_store_rates_peps[sample][cut_name][pipeline] = val_rate_tryptic_pep
                    samples_store_rates_kmers[sample][cut_name][pipeline] = val_rate_kmers
                    
                    print(f'{len( samples_store_kmers[sample][cut_name][pipeline])} validated kmers')

                    print('\n')
                    
    # Performs sets comparisons
    compare_kmers = compare_OHSU_ETH(samples_store_kmers, read_from_disk=True)
    compare_peptides = compare_OHSU_ETH(samples_store_pep, read_from_disk=True)
    peptide_rates = format_validation_rates(samples_store_rates_peps, read_from_disk=True)
    kmers_rates = format_validation_rates(samples_store_rates_kmers, read_from_disk=True)
    
    path_data_pep = os.path.join(save_folder, f'data_peptides{MS_FDR}_{MS_strategy}.tsv.gz')
    path_data_kmers = os.path.join(save_folder, f'data_kmers{MS_FDR}_{MS_strategy}.tsv.gz')
    path_data_peptide_rates = os.path.join(save_folder, f'data_peptides-rates{MS_FDR}_{MS_strategy}.tsv.gz')
    path_data_kmers_rates = os.path.join(save_folder, f'data_kmers-rates{MS_FDR}_{MS_strategy}.tsv.gz')

    compare_peptides.to_csv(path_data_pep, sep = '\t', index = None)
    print(f'Saved data to {path_data_pep}')
    
    compare_kmers.to_csv(path_data_kmers, sep = '\t', index = None)
    print(f'Saved data to {path_data_kmers}')
    
    peptide_rates.to_csv(path_data_peptide_rates, sep = '\t', index = None)
    print(f'Saved data to {path_data_peptide_rates}')
    
    kmers_rates.to_csv(path_data_kmers_rates, sep = '\t', index = None)
    print(f'Saved data to {path_data_kmers_rates}')
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Loads proteomics results and extract validated kmers and peptides for plotting')
    parser.add_argument("--proteomicsdir", help='path with proteomics result (OSHU and ETH) tree')
    
    parser.add_argument("--samples-breast",help='list of breast sample names', nargs='+')
    parser.add_argument("--samples-ov", help='list of breast sample names', nargs='+')
    
    parser.add_argument("--fasta-base-OHSU", help='folder containing fasta files for OHSU')
    parser.add_argument("--fasta-base-ETH", help='folder containing fasta files for ETH')
    
    parser.add_argument("--kmer-files-OHSU", help='tar archiv file containing filtered results files for OHSU')
#     parser.add_argument("--basedir-breast", help='path with filtering result tree for ETH breast samples')
#     parser.add_argument("--basedir-ov", help='path with filtering result tree for ETH ov samples')

    parser.add_argument("--pipelines", help='name of the pipelines in the proteomics result tree', nargs='+', default=['ETH', 'OHSU'])
    
    parser.add_argument("--FDR-limit", type=float, help='FDR limit', default=0.05)
    parser.add_argument("--MS-FDR", help='FDR method: Either peptide FDR with crema or psm FDR with crux', choices=['_crema', '_crux', ''])
    parser.add_argument("--MS-strategy", help='FDR setup: Either single pipeline, or pool per pipeline, or joint pipeline and experiments', choices=['pool', 'joint', 'single'])

    parser.add_argument("--save-folder", help='base folder to save results')
    args = parser.parse_args()
    print(args)
    
    process_proteomics_results(args.proteomicsdir, 
                               args.samples_breast, args.samples_ov, 
#                               args.basedir_breast, args.basedir_ov, 
                               args.fasta_base_OHSU, args.fasta_base_ETH, args.kmer_files_OHSU, 
                               args.pipelines, 
                               args.MS_FDR, args.MS_strategy, args.FDR_limit, 
#                               args.filtering_path_suffix, 
                               args.save_folder)
