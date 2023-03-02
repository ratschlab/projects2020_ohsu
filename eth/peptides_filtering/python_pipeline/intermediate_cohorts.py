#!/usr/bin/env python
# coding: utf-8

import argparse
import numpy as np
import os
import timeit
import glob 
import pandas as pd
import time
import multiprocessing as mp 
import logging
import sys 
from pathlib import Path
from helpers_filter import process_libsize
from helpers_filter import normalization


##### MAIN #####
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run specifications')
    parser.add_argument('--base-cancer', type=str, required=True, help='folder with cancer results')
    parser.add_argument('--base-normal', type=str, required=True, help='folder with normal results')
    parser.add_argument('--cancer-targets', nargs='+', type=str, required=False, default=None, help='Name of target cancer samples')
    parser.add_argument('--interm-cancer-cohort', type=str, required=True, help='basename of intermediate cancer files to merge')
    parser.add_argument('--interm-normal-cohort', type=str, required=True, help='basename of intermediate normal files to merge')
    parser.add_argument('--path-cancer-libsize', type=str, required=True, help='file with cancer libsize')
    parser.add_argument('--normalizer-libsize', type=int, required=True, help='normalizer')
    parser.add_argument('--metadata', nargs='+', type=str, required=False, default=None, help='Metadata of intermediate files')
    parser.add_argument('--intermediate-output', type=str, required=True, help='name of output file with merged ready to filter matrix'
    args = parser.parse_args()

    expr_matrix = 'ref_graph_kmer_JuncExpr'
    sample_lim = 1
    base_cancer = args.base-cancer
    base_normal = args.base-normal 
    cancer_targets = args.cancer-targets
    interm_cancer_cohort = args.interm-cancer-cohort
    interm_normal_cohort = args.interm-normal-cohort
    path_libsize = args.path-cancer-libsize
    normalizer_libsize = args.normalizer-libsize
    metadata = args.metadata
    intermediate_output = args.intermediate-output

    # Cancer cohort files - 3 min 
    start_time  = timeit.default_timer()
    cohort_cancer = glob.glob(os.path.join(base_cancer, 'cohort_mutNone/tmp_out_ref_batch_*', interm_cancer_cohort))
    time_res = timeit.default_timer() - start_time 
    print(f'{len(cohort_cancer)} cohort cancer file')
    print(time_res)

    # Cancer all raw files 
    start_time  = timeit.default_timer()
    path_cohort = glob.glob(os.path.join(base_cancer, 'cohort_mutNone/tmp_out_ref_batch_*'))
    time_res = timeit.default_timer() - start_time 
    print(time_res)


    # GTEX cohort files 4.5 min
    start_time  = timeit.default_timer()
    cohort_gtex = glob.glob(os.path.join(base_gtex,'cohort_mutNone/tmp_out_ref_batch_*', interm_gtex_cohort)) #path_gtex_cohort
    time_res = timeit.default_timer() - start_time 
    print(f'{len(cohort_gtex)} cohort cancer file')
    print(time_res)


    # Annot
    start_time  = timeit.default_timer()
    annot_cancer = glob.glob(os.path.join(base_cancer, 'cohort_mutNone/tmp_out_ref_batch_*/ref_annot_kmer.gz'))
    time_res = timeit.default_timer() - start_time 
    print(f'{len(annot_cancer)} cohort cancer file')
    print(time_res)


    # Out directory
    outdir = os.path.join(base_cancer, 'filtered_cancer') 
    Path(outdir).mkdir(parents=True, exist_ok=True)
    print(f'Creating directory {outdir}')


    # Sample output dir
    for target_sample in target_samples:
        outdir_sample = os.path.join(outdir, target_sample)
        Path(outdir).mkdir(parents=True, exist_ok=True)
        print(f'Creating directory {outdir_sample}')


    # # Create large intermediate table

    # ## Load

    df_expr = []
    report_count = [] 
    report_steps = []


    # Load libsize - OK
    libsize = process_libsize(path_libsize, normalizer_libsize) 


    # Load cancer sample - 20 min 

    start_time  = timeit.default_timer()
    for idx, batch_gene in enumerate(path_cohort):
        if os.path.exists(os.path.join(batch_gene , expr_matrix)) and \
           os.path.exists(os.path.join(batch_gene , 'output_sample_IS_SUCCESS')):
            
            partitions = glob.glob(os.path.join(batch_gene , expr_matrix) + '/*')
            #print(os.path.join(batch_gene , expr_matrix))
            #print(f'Number of partitions is {len(partitions)}')
            for part in partitions:
                kmers_sample = pd.read_csv(part, sep = '\t', usecols = target_samples + metadata)
                kmers_sample['batch'] = os.path.basename(batch_gene).split('_')[-1]
                df_expr.append(kmers_sample)
            #print(f'processed gene {idx}')


    cancer_targets = pd.concat(df_expr, axis = 0)

    time_res = timeit.default_timer() - start_time 
    print(time_res)


    metadata = metadata + ['batch']


    # Normalize sample of interest 
    sample_cols = target_samples 
    cancer_targets = normalization(cancer_targets, sample_cols, libsize, metadata )


    # Set of target kmers
    target_kmers = set(cancer_targets['kmer'])


    print(f'Unique cancer kmers { len(target_kmers)}')


    print(f'Loading: Size cancer targets: { cancer_targets.shape}')


    print(cancer_targets.head())


    # restrict expressed in target samples - hardcoded number of targets

    cancer_targets = cancer_targets.loc[(cancer_targets[target_samples[0]] > 0) | \
                      (cancer_targets[target_samples[1]] > 0) | \
                      (cancer_targets[target_samples[2]] > 0) |\
                      (cancer_targets[target_samples[3]] > 0) | \
                      (cancer_targets[target_samples[4]] > 0) ]


    print(f'Restricting to 5 samples: {cancer_targets.shape}')


    cancer_targets = cancer_targets.drop_duplicates() # drop duplicates



    print(f'Drop Duplicates: cancer target size {cancer_targets.shape}')


    # Load annotation - OK 12 min 
    # Create a set of annotated kmers
    print('Processing annotation')
    start_time  = timeit.default_timer()
    cancer_targets_annotated = set()
    for idx, path in enumerate(annot_cancer):
        annot = pd.read_csv(path, sep = ',') 
        annot = set(annot['kmer'])
        cancer_targets_annotated.update(target_kmers.intersection(annot))
    time_res = timeit.default_timer() - start_time 
    print(time_res)


    print(f'Loading: Annotated unique kmers: {len(cancer_targets_annotated)}')


    # Load cancer cohort - 11 minutes
    start_time  = timeit.default_timer()
    df_cancer_cohort = []
    print(f'Reading cohort cancer')
    for path in cohort_cancer:
        tmp_cancer = pd.read_csv(path, sep = ',') 
        tmp_cancer = tmp_cancer.loc[tmp_cancer['isCrossJunction'] == True]
        tmp_cancer['batch'] = path.split('/')[-2].split('_')[-1]
        df_cancer_cohort.append(tmp_cancer)
    time_res = timeit.default_timer() - start_time 
    print(time_res)


    print(f'Loading: Cancer cohort intermediate file size {len(df_cancer_cohort)}')


    # Load GTEX cohort - OK 1h50
    start_time  = timeit.default_timer()
    df_gtex_cohort = []
    print(f'Reading cohort gtex')
    for idx, path in enumerate(cohort_gtex):
        tmp_gtex = pd.read_csv(path, sep = ',')
        tmp_kmer = set(tmp_gtex['kmer'])
        tmp_kmer = tmp_kmer.intersection(target_kmers)
        tmp_gtex = tmp_gtex.set_index('kmer').loc[tmp_kmer].reset_index()
        df_gtex_cohort.append(tmp_gtex)
    time_res = timeit.default_timer() - start_time 
    print(time_res)


    print(f'Loading: Normal cohort intermediate file size {len(df_gtex_cohort)}')


    # ## Merge

    print(f'Prior to merge: cancer targets size {cancer_targets.shape}')

    df_gtex_cohort = pd.concat(df_gtex_cohort, axis = 0)

    df_cancer_cohort = pd.concat(df_cancer_cohort, axis = 0)

    df_gtex_cohort = df_gtex_cohort.drop_duplicates()

    df_cancer_cohort = df_cancer_cohort.drop_duplicates()


    print(f'Concat and drop duplicates: Normal cohort {len(df_gtex_cohort)}')


    print(f'Concat and drop duplicates: Cancer cohort {len(df_cancer_cohort)}')


    # cohorts
    df_gtex_cohort = df_gtex_cohort.rename({col: 'gtexCohort' + col for col in df_gtex_cohort.columns if 'filter' in col}, axis = 1)

    df_cancer_cohort = df_cancer_cohort.rename({col: 'cancerCohort' + col for col in df_cancer_cohort.columns if 'filter' in col}, axis = 1)


    # annotation
    cancer_targets_annotated = pd.DataFrame(cancer_targets_annotated, columns = ['kmer'])

    cancer_targets_annotated['isAnnotated'] = 1
    cancer_targets_annotated = cancer_targets_annotated.drop_duplicates()


    # Merge annotation 
    cancer_targets = cancer_targets.merge(cancer_targets_annotated, on = 'kmer', how = 'left')


    print(f'After merge annotation: cancer targets size {cancer_targets.shape}')


    # Merge cancer on coord and kmer col
    cancer_targets = df_cancer_cohort.loc[:, metadata + \
                         [col for col in df_cancer_cohort.columns \
                          if 'filter' in col]].merge(cancer_targets , on = metadata, how = 'right')


    print(f'After merge Cancer cohort: cancer targets size {cancer_targets.shape}, merge on metadata columns')


    print(cancer_targets.head())


    # Merge normals on kmer col
    cancer_targets = df_gtex_cohort.loc[:, ['kmer'] + \
                         [col for col in df_gtex_cohort.columns \
                          if 'filter' in col]].merge(cancer_targets , on = ['kmer'] , how = 'right')


    print(f'After merge Normal cohort: cancer targets size {cancer_targets.shape}, merge on kmer columns')

    print(cancer_targets.head())



    print(f'Save intermediate table to {intermediate_output}')

    cancer_targets.to_csv(intermediate_output, compression='gzip', sep = '\t', index=None)

