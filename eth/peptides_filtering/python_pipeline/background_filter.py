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
    whitelist_normal_tag = '10-21overlap'
    whitelist_normal = '/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist'#", help="file containg whitelist for normal samples", required=False, default=None)
    path_normal_libsize = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/ARCHIV_keep_runs/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/expression_counts.libsize.tsv' #help="libsize file path for normal samples", required=False, default=None)
    normalizer_normal_libsize = 400000
    filters = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0] 
    sample_pattern = 'SRR'
    metadata = ['kmer', 'coord', 'junctionAnnotated', 'readFrameAnnotated', 'isCrossJunction']
    # ---------------------------

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
    
    expr_matrix = 'ref_graph_kmer_SegmExpr'
    if os.path.exists(os.path.join(batch_gene , expr_matrix)) and \
       os.path.exists(os.path.join(batch_gene , 'output_sample_IS_SUCCESS')):
        df_gene_batch_filt_Segm, n_partitions = filter_on_partition(os.path.join(batch_gene , expr_matrix),
                                                                   n_partitions, 
                                                                   libsize, 
                                                                   whitelist, 
                                                                   sample_pattern, 
                                                                   metadata,
                                                                   filters)
    else: 
        df_gene_batch_filt_Segm = None
    expr_matrix = 'ref_graph_kmer_JuncExpr'
    if os.path.exists(os.path.join(batch_gene , expr_matrix)) and \
       os.path.exists(os.path.join(batch_gene , 'output_sample_IS_SUCCESS')):
        df_gene_batch_filt_Junc, n_partitions = filter_on_partition(os.path.join(batch_gene , expr_matrix),
                                                                   n_partitions, 
                                                                   libsize, 
                                                                   whitelist, 
                                                                   sample_pattern, 
                                                                   metadata,
                                                                   filters)
    else:
        df_gene_batch_filt_Junc = None
        
    if (df_gene_batch_filt_Segm is not None) and (df_gene_batch_filt_Junc is not None): 
        res = pd.concat([df_gene_batch_filt_Segm, df_gene_batch_filt_Junc], axis = 0)
    elif (df_gene_batch_filt_Segm is not None) and (df_gene_batch_filt_Junc is None):
        res = df_gene_batch_filt_Segm
    elif (df_gene_batch_filt_Segm is None) and (df_gene_batch_filt_Junc is not None):
        res = df_gene_batch_filt_Junc
    else:
        res = None

    
    if res is not None:
        outfile = os.path.join(batch_gene, f'ref_graph_kmer{tag_normalize}filtered{whitelist_normal_tag}.gz')
        res.to_csv(outfile, compression = 'gzip', index = None)
        #now = datetime.now()
        #current_time = now.strftime("%H:%M:%S")
        current_time=''
        print(f'{current_time}: Saved to {outfile}', flush=True)
    return 'done'

def handler(error):
    print(f'Error: {error}', flush=True)
    
def dummy(path_cohort):
    print(path_cohort, flush=True)
    try:
        print('hello')
    except:
        print('NONE')

##### MAIN #####
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run specifications')
    parser.add_argument('--processes', type=int, required=True, help='the number of processes for the multiprocessing')
    args = parser.parse_args()
    
    path_cohort = glob.glob('/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/GTEX2019_eth/GTEX2019_c4dd02c_conf2_RFall_ref/cohort_mutNone/*')
    path_cohort=path_cohort[920:1000]
    #path_cohort = path_cohort[0:10] #TODO remove
    print(f'{len(path_cohort)} batches found', flush=True)
    print(f'Run with {args.processes} processes')

    pool = mp.Pool(args.processes)
    result = pool.map_async(process_on_cohort, path_cohort,  error_callback=handler, chunksize=2) 
    result.wait()
    pool.close()
    pool.join()


