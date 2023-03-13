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



def process_on_cohort(batch_gene, whitelist, whitelist_tag, path_libsize, normalizer, filters, sample_pattern, do_overwrite=False, do_normalize=True):
    ''' Main code for filtering'''
    
    metadata = ['kmer', 'coord', 'junctionAnnotated', 'readFrameAnnotated', 'isCrossJunction']
    # ---------------------------
    n_partitions = 0 #TODO make across processes
     
    # Process library size file
    if do_normalize:
        tag_normalize = '_normalized_'
        libsize = process_libsize(path_libsize, normalizer) 

    else:
        tag_normalize = ''
        libsize = None
    
    # Process whitelist
    whitelist, whitelist_tag = process_whitelist(whitelist, whitelist_tag) 
    
    # Define outputs
    outfile = os.path.join(batch_gene, f'ref_graph_kmer{tag_normalize}filtered{whitelist_tag}.gz')

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
    parser.add_argument('--path-cohort', type=str, required=False, default=None, help='immunopepper output path containing the batch folders')
    parser.add_argument('--whitelist-tag', type=str, required=False, default=None, help='suffix to add to the output file') 
    parser.add_argument('--whitelist', type=str, required=False, default=None, help='path to the whitelist of samples. File without header with sample names')
    parser.add_argument('--path-libsize', type=str, required=False, default=None, help='path to the library size. First column sample, Second column libsize_75percent') 
    parser.add_argument('--normalizer-libsize', type=float, required=False, default=None, help='factor by which to multiply after the division by the library size. Default median of libsize')
    parser.add_argument('--filters', nargs='+', type=float, required=False, default=None, help='Thresholds for filtering (list of floats). The number of samples >= (or >0) to each threshold will be computed')
    parser.add_argument('--sample-pattern', type=str, required=False, default=None, help='Commun naming pattern shared by all the sample names') 
    parser.add_argument('--do-overwrite', action="store_true", required=False, default=False, help='Overwrite the outputs. Otherwise, missing outputs in the corresponding batches')
    parser.add_argument('--do-normalize', action="store_true", required=False, default=False, help='If set, perform a normalisation with the library size')
    args = parser.parse_args()
    
    path_cohort = glob.glob(args.path_cohort + '/*')
 #   whitelist_tag = '10-21overlap_TEST'
 #   whitelist = '/cluster/work/grlab/projects/projects2020_OHSU/sample_lists/GTEX/GTEx_sample_IDs_10-2021_lib_graph_juliannelist'#", help="file containg whitelist for normal samples", required=False, default=None)
 #   path_libsize = '/cluster/work/grlab/projects/TCGA/PanCanAtlas/immunopepper_paper/peptides_ccell_rerun_gtex_151220/ARCHIV_keep_runs/GTEX2019_commit_v3_TEST_merged3_372a147_medium_run_pya.0.17.1_conf2_annot_ref_chrall_cap/expression_counts.libsize.tsv' #help="libsize file path for normal samples", required=False, default=None)
 #   normalizer_libsize = 400000
 #   filters = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0]
 #   sample_pattern = 'SRR'
 #   do_overwrite = False  # do not overwrite computed files
 #   do_normalize = True


    if args.start_id is not None and args.end_id is not None: 
        path_cohort=path_cohort[args.start_id:args.end_id]
    elif args.start_id is not None:
        path_cohort=path_cohort[args.start_id:]
    elif args.end_id is not None:
        path_cohort=path_cohort[:args.end_id]
    input_args = [(batch_gene, args.whitelist, args.whitelist_tag, args.path_libsize, args.normalizer_libsize, args.filters, args.sample_pattern, args.do_overwrite, args.do_normalize) for batch_gene in path_cohort]
    print(args, flush=True)
    print(f'Start id: {args.start_id} and End id: {args.end_id}', flush=True)
    print(f'{len(path_cohort)} batches found', flush=True)
    print(f'Run with {args.processes} processes', flush=True)
    print(f'---- Starting multiprocessing ----', flush=True)
    pool = mp.Pool(args.processes)
    result = pool.starmap_async(process_on_cohort, input_args, error_callback=handler, chunksize=2) 
    result.wait()
    pool.close()
    pool.join()


