import pandas as pd
import os
import glob
import timeit
import numpy as np
import h5py
from collections import defaultdict



def explode_immunopepper_coord(mx):
    coord_mx = mx['coord'].str.split(':', expand=True) #7 min

    coord_mx[1] = coord_mx[1].astype('int')
    coord_mx[2] = coord_mx[2].astype('int')

    coord_mx['strand'] = None
    coord_mx.loc[coord_mx[1] < coord_mx[2], 'strand'] = '+'
    coord_mx.loc[coord_mx[1] > coord_mx[2], 'strand'] = '-'

    coord_mx['junction_coordinate1'] = None
    coord_mx['junction_coordinate2'] = None


    coord_mx = coord_mx.astype(str) # 7 min

    coord_mx['+first'] = coord_mx[1] + ':' + coord_mx[2]
    coord_mx['+secon'] = coord_mx[3] + ':' + coord_mx[4]
    coord_mx['-first'] = coord_mx[3] + ':' + coord_mx[0]
    coord_mx['-secon'] = coord_mx[5] + ':' + coord_mx[2]

    coord_mx.loc[(coord_mx['strand'] == '+'), 'junction_coordinate1'] = coord_mx['+first'] 
    coord_mx.loc[(coord_mx['strand'] == '-'), 'junction_coordinate1'] = coord_mx['-first'] 
    coord_mx.loc[(coord_mx['strand'] == '+') & (coord_mx[4] != 'None') , 'junction_coordinate2'] = coord_mx['+secon']
    coord_mx.loc[(coord_mx['strand'] == '-') & (coord_mx[4] != 'None') , 'junction_coordinate2'] = coord_mx['-secon']
    
    return coord_mx


def read_libsize_whitelist(libsize, whitelist):
    # Read libsize and whitelist
    libsize_normal = pd.read_csv(libsize, sep = '\t')
    whitelist_normal = pd.read_csv(whitelist, sep = '\t', header = None)
    whitelist_normal.columns = ['sample']
    return libsize_normal, whitelist_normal


def preprocess_STAR_projected(chrm, path_star, whitelist_normal, libsize_normal):
    # Star junctions - projected coordinates and expression
    projected_chr_file = os.path.join(path_star, f'junctions_spladder_projected/junctions_spladder.projected.{chrm}.hdf5')
    star_expr = h5py.File(projected_chr_file, 'r')

    # Whitelist samples
    index_whitelist_samples = [s for s, sample in enumerate(star_expr['samples']) if sample.decode() + 'all' in whitelist_normal['sample'].values]
    samples_decoded = [sample.decode() + 'all' for s, sample in enumerate(star_expr['samples']) if sample.decode() + 'all' in whitelist_normal['sample'].values]

    # Libsize sampple     
    lib_75_per_sample = libsize_normal.set_index('sample').loc[samples_decoded, 'libsize_75percent'].values

    return star_expr, index_whitelist_samples, lib_75_per_sample



def h5_indexes(junction_start, junction_end, chrm, strand, expression_h5):
    # Extract junction ID in hdf5 file 
    jx_idx_h5 = np.where((junction_start == expression_h5[f'{chrm}:{strand}:junction_start'][...]) & 
            (junction_end == expression_h5[f'{chrm}:{strand}:junction_end'][...]))[0]
    return jx_idx_h5


def h5_expression(expression_h5, chrm, strand, jx_idx_h5, index_whitelist_samples, 
                  lib_75_per_sample, normalizer):
    # Extract Raw count in hdf5 file
    #start_time = timeit.default_timer()    
    raw_counts = expression_h5[f'{chrm}:{strand}:count'][jx_idx_h5, :]
    raw_counts = raw_counts[np.array(index_whitelist_samples)]
    #print(timeit.default_timer() - start_time, 'raw counts')

    normalized_counts = np.divide(raw_counts, lib_75_per_sample) * normalizer
    return normalized_counts


def get_junction_counts(junction_start, junction_end, 
                        chrm, strand, expression_h5, 
                        index_whitelist_samples, lib_75_per_sample, normalizer):
    strand_complementary = {'+':'-', '-':'+'}
    
    # Indexes 
    jx_idx_h5 = h5_indexes(junction_start, junction_end, chrm, strand, expression_h5)
    if len(jx_idx_h5) == 0:
        strand = strand_complementary[strand]
        jx_idx_h5 = h5_indexes(junction_start, junction_end, chrm, strand, expression_h5)

    # Expression     
    if len(jx_idx_h5 >0):
        assert(len(jx_idx_h5) == 1) #TODO check critical
        jx_idx_h5 = jx_idx_h5[0]
        normalized_counts = h5_expression(expression_h5, chrm, strand, 
                                          jx_idx_h5, index_whitelist_samples, lib_75_per_sample, normalizer)
    else:
        normalized_counts = None
        print(f'Error: No {chrm}:strand:junction_start {chrm}:strand:junction_end matching {junction_start}:{junction_end}', flush=True)
    return normalized_counts


def filter_recurrence(array_, threshold):
    if threshold == 0:
        return np.sum(array_ > threshold)
    else:
        return np.sum(array_ >= threshold)
    
    
def filter_multi_thresholds(normalized_counts, filter_thresholds):
    recurrence = []
    for threshold in filter_thresholds: # Make filter threshold
        recurrence.append(filter_recurrence(normalized_counts, threshold))
    return recurrence 


def collect_expression_thresholds(libsize, whitelist, normalizer, filter_thresholds, 
                                  path_star, chr_jx, jx_strand ):
    res = []
    libsize_normal, whitelist_normal = read_libsize_whitelist(libsize, whitelist)
    counter = 0
    for chrm, jxS in chr_jx.items(): # Per chromosome
        print(f'{chrm}: with {len(jxS)} junctions', flush=True)
        # Expression file
        expression_h5, index_whitelist_samples, lib_75_per_sample = preprocess_STAR_projected(chrm, 
                                                                                          path_star, 
                                                                                          whitelist_normal, 
                                                                                          libsize_normal)
        start_time = timeit.default_timer()
        for jx in jxS: # Per junction
            counter +=1
            if counter % 500 == 0:
                print(f'....{counter}', flush=True)
            junction_start = int(jx.split(':')[0])
            junction_end = int(jx.split(':')[1])
            normalized_counts = get_junction_counts(junction_start, junction_end, 
                                                    chrm, jx_strand[jx], expression_h5, 
                                                    index_whitelist_samples, 
                                                    lib_75_per_sample, normalizer)

            if normalized_counts is not None:
                metadata = [jx, jx_strand[jx], chrm] + filter_multi_thresholds(normalized_counts, filter_thresholds)

            res.append(metadata)

        expression_h5.close() 
        print(f'{timeit.default_timer() - start_time} seconds', flush=True)
        
    return res
