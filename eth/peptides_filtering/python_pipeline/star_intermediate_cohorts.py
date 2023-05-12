#!/usr/bin/env python
# coding: utf-8



import pandas as pd
import os
import glob
import timeit
import numpy as np
import h5py
from collections import defaultdict
from helpers_STAR import * 
import argparse


##### MAIN #####
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='run specifications')
    parser.add_argument('--path-star', type=str, required=True, help='folder with star junctions')
    parser.add_argument('--big-matrix', type=str, required=True, default=None, help='Intermediate filtering matrix')
    parser.add_argument('--whitelist', type=str, required=True, help='whitelist of samples to include in recurrence analysis')
    parser.add_argument('--libsize', type=str, required=True, help='library size file')
    parser.add_argument('--normalizer', type=int, required=True, help='value to rescale after normalisation')
    parser.add_argument('--filter-thresholds', nargs='+', type=float, required=True, default=None, 
                        help='Number of reads to use as expression threshold in recurrence analysis')
    parser.add_argument('--jx-target-list', type=str, required=False, default='', help='path to file with whitelist kmers coordinates')
    args = parser.parse_args()

    print(args, flush=True)


    # -------- Paths -------- 

    # -------- Path STAR
    path_star = args.path_star
    # Star junctions - unique coordinates
    star_jx = os.path.join(path_star, 'junctions_spladder.all_coords.sorted.uniq.tsv.gz')
    # Star junctions - projected coordinates and expression
    #projected_chr_file = os.path.join(path_star, f'junctions_spladder_projected/junctions_spladder.projected.{chrm}.hdf5')


    # -------- Intermediate filtering results (threshold and merged)
    #TODO DO the generation matrix???
    # Foreground matrix
    big_matrix = args.big_matrix



    # -------- GTEX filtering 
    whitelist = args.whitelist
    libsize = args.libsize


    normalizer = args.normalizer
    filter_thresholds = args.filter_thresholds
    #Optional if needs to process a junction list
    jx_target_list = args.jx_target_list




    # -------- Load -------- 

    # Load cancer generation matrix
    mx = pd.read_csv(big_matrix, sep = '\t')
    #display(mx.head())
    # Add split junction information to generation table
    coord_mx = explode_immunopepper_coord(mx)
    #display(coord_mx.head())
    mx = pd.concat([mx, coord_mx[['strand', 'junction_coordinate1', 'junction_coordinate2']]], axis = 1)
    print('foreground matrix shape', mx.shape, flush=True)

    if jx_target_list:
        target_jx = pd.read_csv(jx_target_list, header = None) 
        target_jx.columns = ['coord']
        mx = mx.merge(target_jx, on = 'coord', how = 'inner')
    print('foreground matrix shape', mx.shape, flush=True)   
        
    # LOAD STAR junctions
    star_jx = pd.read_csv(star_jx, sep = '\t')
    star_jx.head()
    # Add STAR junction column 10 min
    star_jx['junction_coordinate'] = star_jx['junction_start'].astype(str) + ':' + star_jx['junction_end'].astype(str)


    # --------  Check junction presence in STAR --------  
    #- 2-exons, 1st junction is present 
    #- 3 exons, needs first and second junction present 

    # Foreground Kmers from 1 junction
    one_jx = mx[(mx['junction_coordinate1'] != 'None') & (mx['junction_coordinate2'] == 'None')]
    print('foreground matrix one junction', one_jx.shape, flush=True)

    # Foreground Kmers from 2 junctions
    two_jx = mx[(mx['junction_coordinate1'] != 'None') & (mx['junction_coordinate2'] != 'None')]
    print('foreground matrix two junctions', two_jx.shape, flush=True)

    # Foreground 1 junction - NOT IN STAR 
    isstar_one = set(one_jx['junction_coordinate1']).intersection(set(star_jx['junction_coordinate'])) # junction coordinates
    print('foreground matrix one junction is in star', len(isstar_one), flush=True)
    is_star_one_jx = set(one_jx.set_index('junction_coordinate1').loc[isstar_one, 'coord']) # corresonding kmer coordinates
    print('-> foreground matrix one-junction-kmers are in star', len(is_star_one_jx), flush=True)

    # Foreground 2 junctions - NOT IN STAR 
    is_star_two1 = set(two_jx['junction_coordinate1']).intersection(set(star_jx['junction_coordinate'])) # Junction coordinates
    is_star_two2 = set(two_jx['junction_coordinate2']).intersection(set(star_jx['junction_coordinate'])) # Junction coordinates
    print('foreground matrix two junctions: left junction is in star', len(is_star_two1), flush=True)
    print('foreground matrix two junctions: right junction is in star', len(is_star_two2), flush=True)

    is_star_two_jx = set(two_jx.set_index('junction_coordinate1').loc[is_star_two1, 'coord']).intersection(\
    set(two_jx.set_index('junction_coordinate2').loc[is_star_two2, 'coord'])) # Corresponding kmer coord
    print('-> foreground matrix two-junctions-kmers are in star', len(is_star_two_jx), flush=True)


    # Foreground Table - Create FLAG GTEX junctions
    start_time = timeit.default_timer()
    mx['STAR_GTEx_jx'] = False
    mx = mx.set_index('coord')
    mx.loc[list(is_star_one_jx.union(is_star_two_jx)), 'STAR_GTEx_jx'] = True
    mx = mx.reset_index()
    print('time create STAR junction presence flag', timeit.default_timer() - start_time, flush=True)    

    # --------  Retrieve expression -------- 

    # Extract metadata

    # Select junctions to query in STAR file 
    query = mx.loc[mx['STAR_GTEx_jx'] == True, ['strand', 'junction_coordinate1',
                                                'junction_coordinate2', 'STAR_GTEx_jx']].drop_duplicates()\
                                                                                        .reset_index()\
                                                                                        .drop('index', axis = 1)

    assert(query.loc[:, ['junction_coordinate1', 'strand']].drop_duplicates().shape[0] == \
           query.loc[:, ['junction_coordinate1']].drop_duplicates().shape[0])

    # Query-Junction to strand 
    start_time = timeit.default_timer()
    jx_strand = {}
    for i in np.arange(query.shape[0]):
        jx_strand[query.loc[i, 'junction_coordinate1']] = query.loc[i, 'strand']

    for i in np.arange(query.shape[0]):
        jx_strand[query.loc[i, 'junction_coordinate2']] = query.loc[i, 'strand']

    jx_strand.pop('None')
    print('Number of query junctions:', len(jx_strand), flush=True)
    print('time Query-junction to strand', timeit.default_timer() - start_time, flush=True)

    # Chromosome to query-junction
    start_time = timeit.default_timer()
    star_jx = star_jx.set_index('junction_coordinate') #Faster
    intermediate = star_jx.loc[jx_strand.keys(), :]

    chr_jx = defaultdict(set)
    for jx, chrm in zip(intermediate.index, intermediate['chr']):
        chr_jx[chrm].add(jx)
    print('chromosomes involved', chr_jx.keys(), flush=True)
    print('Chromosome to query-junction', timeit.default_timer() - start_time, flush=True)
        

    ### Collect Expression 

    # -------- Thresholding experiments -------- 
    # Remark 1 junction can have multiple chrm --> Which 
    # Remark 1 junction on 1 chrm can have multiple strands --> Take from Immunopepper

    print('Start collect expression', flush=True)
    res = collect_expression_thresholds(libsize, whitelist, normalizer,
                                        filter_thresholds, path_star, chr_jx, jx_strand )



    gtex_cols = ['gtexCohortfilter >0.0', 'gtexCohortfilter >=1.0',
           'gtexCohortfilter >=2.0', 'gtexCohortfilter >=3.0',
           'gtexCohortfilter >=5.0', 'gtexCohortfilter >=10.0']

    # Expression to DF
    df_res = pd.DataFrame(res, columns = ['junction_coordinate', 'strand_STAR', 'chr'] + gtex_cols)
    print('Junction - expressin df shape', df_res.shape, flush=True)
    #display(df_res.head())

    #-------- Format results -------- 



    # One junction merge (straightforward merge on junction)
    print('foreground matrix one junction - BIS', one_jx.shape, flush=True)
    # Remove the GTEX quantifications from immunopepper
    base_one_jx = one_jx.drop(gtex_cols, axis = 1).drop_duplicates() 
    base_one_jx.head()
    print('foreground matrix one junction no GTEx', base_one_jx.shape, flush=True)
    # Add GTEX quantifications from STAR
    one_jx_quantified = base_one_jx.merge(df_res, left_on = 'junction_coordinate1', 
                                          right_on = 'junction_coordinate', how = 'left') 
    print(one_jx_quantified.columns, flush=True)
    #display(one_jx_quantified.tail())

    # Two junctions merge (not straight forwrds, consider botrh junctions separately and take the max)
    print('foreground matrix two junctions - BIS', two_jx.shape, flush=True)
    # Remove the GTEX quantifications from immunopepper
    base_two_jx = two_jx.drop(gtex_cols, axis = 1).drop_duplicates() 
    base_two_jx.head()
    print('foreground matrix two junctions no GTEx', base_two_jx.shape, flush=True)
    # Add GTEX quantifications from STAR  # first junction
    two_jx_quantified_left = base_two_jx.merge(df_res, left_on = 'junction_coordinate1', 
                                          right_on = 'junction_coordinate', how = 'left')


    # Add GTEX quantifications from STAR  # second junction
    two_jx_quantified_right = base_two_jx.merge(df_res, left_on = 'junction_coordinate2', 
                                          right_on = 'junction_coordinate', how = 'left')

    col_merge = [col for col in two_jx_quantified_right if col not in gtex_cols]

    # Add GTEX quantifications from STAR  # both junctions
    two_jx_quantified = two_jx_quantified_left.merge(two_jx_quantified_right, on = col_merge, how = 'outer')

    # Add GTEX quantifications from STAR  # MAX (!!) over 2 junctions
    for col in gtex_cols:
        two_jx_quantified[col] = two_jx_quantified[[col + '_x', col + '_y']].max(skipna = True, axis = 1)
        two_jx_quantified = two_jx_quantified.drop([col + '_x', col + '_y'], axis = 1)
        
    print(two_jx_quantified.columns, flush=True)
    #display(two_jx_quantified.tail())

    # All Junctions quantified on STAR for GTEX
    print('foreground matrix one junction merged with expression', one_jx_quantified.shape, flush=True)
    print('foreground matrix two junctions merged with expression', two_jx_quantified.shape, flush=True)

    jx_quantified = pd.concat([one_jx_quantified, two_jx_quantified])
    print('foreground matrix all merged with expression', jx_quantified.shape, flush=True)

    ### Save 
    save_path = big_matrix.replace('tsv.gz', 'STAR_GTEx.tsv.gz')
    print('Saving to', save_path, flush=True)
    jx_quantified.to_csv(save_path, compression = 'gzip', index = False, sep = '\t')

