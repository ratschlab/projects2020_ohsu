import argparse
from datetime import datetime
import glob 
from helpers_filter import *
import multiprocessing as mp 
import numpy as np
import os
import pandas as pd
import timeit
import time



def process_on_cohort(batch_gene):
    ''' Main code for background cohort filtering'''
    
    # --- Background Specific ---
    # Will be run once, hence no proper command line 
    whitelist_normal_tag = ''
    whitelist_normal = '/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/TCGA_foreground/sample_full_BRCA_1102_format.tsv'
    path_normal_libsize = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/expression_counts.libsize.tsv'
    normalizer_normal_libsize = 400000
    filters = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0] 
    sample_pattern = 'TCGA'
    metadata = ['kmer', 'coord', 'junctionAnnotated', 'readFrameAnnotated', 'isCrossJunction']
    # ---------------------------
    do_overwrite=True  # do not overwrite computed files
    do_normalize = True
    n_partitions = 0 #TODO make across processes
     
    # Process library size file
    if do_normalize:
        tag_normalize = '_normalized_'
        libsize = process_libsize(path_normal_libsize, normalizer_normal_libsize) #TODO global

    else:
        tag_normalize = ''
        libsize = None
    
    # Process whitelist
    whitelist, whitelist_normal_tag = process_whitelist(whitelist_normal, whitelist_normal_tag) #TODO global
    
    # Define outputs
    outfile = os.path.join(batch_gene, f'ref_graph_kmer{tag_normalize}filtered{whitelist_normal_tag}.gz')

    expr_matrix = 'ref_graph_kmer_SegmExpr'
    df_gene_batch_filt_Segm = None
    df_gene_batch_filt_Junc = None
    if os.path.exists(os.path.join(batch_gene , expr_matrix)) and \
       os.path.exists(os.path.join(batch_gene , 'output_sample_IS_SUCCESS')):
        if do_overwrite or (not os.path.exists(outfile)):
            df_gene_batch_filt_Segm, n_partitions = filter_on_partition(os.path.join(batch_gene , expr_matrix),
                                                                        n_partitions,
                                                                        libsize, 
                                                                        whitelist, 
                                                                        sample_pattern, 
                                                                        metadata,
                                                                        filters)
    expr_matrix = 'ref_graph_kmer_JuncExpr'
    if os.path.exists(os.path.join(batch_gene , expr_matrix)) and \
       os.path.exists(os.path.join(batch_gene , 'output_sample_IS_SUCCESS')):
        if do_overwrite or (not os.path.exists(outfile)):
            df_gene_batch_filt_Junc, n_partitions = filter_on_partition(os.path.join(batch_gene , expr_matrix),
                                                                        n_partitions, 
                                                                        libsize, 
                                                                        whitelist, 
                                                                        sample_pattern, 
                                                                        metadata,
                                                                        filters)
        
    if (df_gene_batch_filt_Segm is not None) and (df_gene_batch_filt_Junc is not None): 
        res = pd.concat([df_gene_batch_filt_Segm, df_gene_batch_filt_Junc], axis = 0)
    elif (df_gene_batch_filt_Segm is not None) and (df_gene_batch_filt_Junc is None):
        res = df_gene_batch_filt_Segm
    elif (df_gene_batch_filt_Segm is None) and (df_gene_batch_filt_Junc is not None):
        res = df_gene_batch_filt_Junc
    else:
        res = None

    
    if res is not None:
        res.to_csv(outfile, compression = 'gzip', index = None)
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        print(f'{current_time}: Saved to {outfile}', flush=True)
    return 'done'

def handler(error):
    print(f'Error: {error}', flush=True)
    

##### MAIN #####
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run specifications')
    parser.add_argument('--processes', type=int, required=True, help='the number of processes for the multiprocessing')
    parser.add_argument('--start-id', type=int, required=False, default=None, help='start id of the subset of batches')
    parser.add_argument('--end-id', type=int, required=False, default=None, help='end id of the subset of batches')
    args = parser.parse_args()
    
    path_cohort = glob.glob('/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102/cohort_mutNone/*')
    if args.start_id is not None and args.end_id is not None: 
        path_cohort=path_cohort[args.start_id:args.end_id]
    elif args.start_id is not None:
        path_cohort=path_cohort[args.start_id:]
    elif args.end_id is not None:
        path_cohort=path_cohort[:args.end_id]

    print(f'Start id: {args.start_id} and End id: {args.end_id}', flush=True)
    print(f'{len(path_cohort)} batches found', flush=True)
    print(f'Run with {args.processes} processes', flush=True)
    print(f'---- Starting multiprocessing ----', flush=True)
    pool = mp.Pool(args.processes)
    result = pool.map_async(process_on_cohort, path_cohort,  error_callback=handler, chunksize=2) 
    result.wait()
    pool.close()
    pool.join()


