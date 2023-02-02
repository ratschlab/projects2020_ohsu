import numpy as np
import os
import timeit
import glob 
import pandas as pd
import time
from IPython.core.debugger import set_trace


def process_libsize(path_lib, custom_normalizer):
    '''
    Loads and returns the normalisation values per sample
    Normalisation formulae: (count / 75 quantile expression in sample) * A
    If no normalisation factor is provided: A = median across samples
    If normalisation factor is provided: A = normalisation factor
    :param path_lib: str path with library size file
    :param custom_normalizer: custum normalisation factor
    :return: dataframe with 75 quantile expression in sample / A values
    '''
    lib = pd.read_csv(path_lib, sep='\t')
    if custom_normalizer:
        lib['libsize_75percent'] = lib['libsize_75percent'] / custom_normalizer
    else:
        lib['libsize_75percent'] = lib['libsize_75percent'] / np.median(lib['libsize_75percent'])
    lib['sample'] = [sample.replace('-', '').replace('.', '').replace('_','') for sample in lib['sample']]
    lib = lib.set_index('sample')
    return lib


def normalization(df, cols, libsize, metadata):
    df_expr = df.loc[:, cols]
    df_norm_vals = libsize.loc[cols, 'libsize_75percent']
    df = pd.concat([df.loc[:, metadata], df_expr / df_norm_vals], axis = 1)
    return df


def process_whitelist(path_whitelist, whitelist_normal_tag):
    if path_whitelist is not None:
        whitelist_samples = set(pd.read_csv(path_whitelist, header = None).iloc[:, 0])
        return whitelist_samples, f'_{whitelist_normal_tag}_'
    else: 
        return None, ''
    
    
def filter_supeq(df, threshold, cols):
    filter_col = f'filter >={threshold}'
    df[filter_col] = np.sum(df[cols] >= threshold, axis = 1)
    return df, filter_col


def filter_supstrict(df, threshold, cols):
    filter_col = f'filter >{threshold}'
    df[filter_col] = np.sum(df[cols] > threshold, axis = 1)
    return df, filter_col



def filter_function(idx, path, libsize, whitelist, sample_pattern, metadata, filters):
    
    ''' 1. Loads the partition 
        2. Normalize, applies the whitelist
        3. Applies all the filters to the partition'''
    
    filter_cols = []
    try:
        start_time  = timeit.default_timer()
        df = pd.read_csv(path, sep = '\t')
        #print(f'Load in {np.round( (timeit.default_timer() - start_time) / 60, 2)} minutes.', flush = True)

        
          
        sample_cols = set([ col for col in df.columns if sample_pattern in col])# --- Background Specific ---
        
        if whitelist:
            sample_cols = sample_cols.intersection(whitelist)   
        sample_cols = list(sample_cols)
        
        start_time  = timeit.default_timer()
        if libsize is not None:
            df = normalization(df, sample_cols, libsize, metadata)
        #print(f'Normalize in {np.round( (timeit.default_timer() - start_time) / 60, 2)} minutes.', flush = True)

        start_time  = timeit.default_timer()
        for read_level in filters:
            if read_level:
                df, col = filter_supeq(df, read_level, sample_cols)
                filter_cols.append(col)
            else: # 0 case
                df, col = filter_supstrict(df, read_level, sample_cols)
                filter_cols.append(col)
        #print(f'Filters in {np.round( (timeit.default_timer() - start_time) / 60, 2)} minutes.', flush = True)

        df = df.loc[:, metadata +  filter_cols]
        return df
    except EOFError: 
         print(f'...Cannot read file {path}. Skipping it.')

            
def filter_on_partition(expr_matrix, n_partitions, libsize, whitelist, sample_pattern, metadata, filters):  
    '''Applies a filtering function to many partitions and concatenate the result'''
    
    # List partitions for gene batch 
    start_time  = timeit.default_timer()
    path_partions = glob.glob(os.path.join(expr_matrix, 'part*'))
    print(f'Processed glob partitions in {np.round( (timeit.default_timer() - start_time) / 60, 2)} minutes.', flush = True)
    N_parts = len(path_partions)
    
    if N_parts:
        n_partitions += N_parts
        start_time  = timeit.default_timer()
        #print(f'\n Start filter: {start_time}, {expr_matrix}', flush = True) 

        df_gene_batch = []
        #print(f'... {N_parts} parts', flush = True)
        for idx, part in enumerate(path_partions):
            start_time  = timeit.default_timer()
            df_gene_batch.append(filter_function(idx, part, libsize, whitelist, sample_pattern, metadata, filters))
            #print(f'Filtered 1 partition in {np.round( (timeit.default_timer() - start_time) / 60, 2)} minutes.', flush = True)
        
        start_time  = timeit.default_timer()
        df_gene_batch = pd.concat(df_gene_batch, axis = 0)   
        #print(f'Concatenation in {np.round( (timeit.default_timer() - start_time) / 60, 2)} minutes.', flush = True)


        time_res = timeit.default_timer() - start_time 
        #print(f'Processed {N_parts} parts in {np.round(time_res/ 60, 2)} minutes. Total partitions so far {n_partitions}', 
        #      flush = True) 
    else: 
        df_gene_batch = None
    
    return df_gene_batch, n_partitions
