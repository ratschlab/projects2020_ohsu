#!/usr/bin/env python
# coding: utf-8

import numpy as np
import os
import timeit
import glob 
import pandas as pd
import time
import multiprocessing as mp 
import logging
import sys 
import pathlib
from pathlib import Path
import matplotlib.pyplot as plt 
from helpers_ffi import * 

# # Order 
# - Uniprot 
# - Sample init
# - Sample expression 
# - Sample cohort 
# - annotation 
# - GTEX?


run_type = 'brca'

# Inputs

if run_type == 'brca':
    target_samples = ['TCGA-AO-A0JM-01A-21R-A056-07.all',
                      'TCGA-C8-A12P-01A-11R-A115-07.all',
                      'TCGA-BH-A18V-01A-11R-A12D-07.all',
                      'TCGA-A2-A0D2-01A-21R-A034-07.all',
                      'TCGA-A2-A0SX-01A-12R-A084-07.all']
    basedir = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Breast_1102'
    intermediate_output = os.path.join(basedir, 'filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz')
elif run_type == 'ov':
    target_samples = ['TCGA-25-1319-01A-01R-1565-13.all',
                      'TCGA-25-1313-01A-01R-1565-13.all',
                      'TCGA-61-2008-01A-02R-1568-13.all',
                      'TCGA-24-1431-01A-01R-1566-13.all',
                      'TCGA-24-2298-01A-01R-1569-13.all']
    basedir = '/cluster/work/grlab/projects/projects2020_OHSU/peptides_generation/CANCER_eth/commit_c4dd02c_conf2_Frame_cap0_runs/TCGA_Ovarian_374'
    intermediate_output = os.path.join(basedir, 'filtering_intermediate/complete_cancer_candidates_order_r.tsv.gz')


# Outputs
filtering_id = 'filters_22March_order_wany_wAnnot'

output_dir = os.path.join(basedir, 'filtering_samples', filtering_id)
pathlib.Path(output_dir).mkdir(exist_ok=True, parents=True)


# Discussion 02/22 Choices
# BACKGROUND cohorts we do (cohort_reads, sample_number)- KEEP pipeline as such
# cohort_reads=[0,1,3, None]
# sample_number=[1,2,10, None]
# FOREGROUND  (cohort_reads, sample_number) means
# cohort_reads=[0,2]
# sample_number(rest of cohort) =[1, 5]
# Sample reads = [0]
# + A case with not foreground


# Parameters
Threshold_target = [0.0]
Threshold_cancer_cohort = [None, 0.0, 2.0, ] # choices = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0]
N_samples_cancer = [None, 1, 5] # choices 1 to 1102 for BRCA and 374 for OV   

Threshold_normal_cohort = [0.0, 1.0, 3.0, None]   # choices = [0.0, 1.0, 2.0, 3.0, 5.0, 10.0]
N_samples_normal = [1, 2, 10, None] #choices 1 to max number of samples in Normal whitelist

tag_cancer = 'cancerCohort'
tag_normal = 'gtexCohort'


tag_prefix = 'G_'
mutation_mode = 'ref'
save_tag = 'GtexCohort'

metadata_save = ['kmer', 'coord', 'junctionAnnotated', 'readFrameAnnotated']

filter_annot = False


# Load matrix to be filtered
df_load = pd.read_csv(intermediate_output, sep = '\t')
print(f'Loaded {intermediate_output}')
df_load = df_load.rename({'batch': f'batch_{run_type}'}, axis = 1)
df_load.shape


df_load.head()
df_load.shape


for cancer_sample_ori in target_samples: # TODO update
    # Sample naming
    target_sample = cancer_sample_ori.replace('-', '').replace('.', '')
    cancer_sample_ori = cancer_sample_ori.replace('.all', '')
    print(f'-------- processing {target_sample} -------- \n')
    
    # Summary file for sample
    summary_file = f'{tag_prefix}filtered_df_{cancer_sample_ori}_samp_chrt_norm_mot.tsv'
    summary_path = os.path.join(output_dir, summary_file)
    print(f'Saving to summary file {summary_path}')
    

    df_expr = []
    report_count = [] 
    report_steps = []
    for threshold_target in Threshold_target:
        for threshold_cancer_cohort in Threshold_cancer_cohort:
            for n_samples_cancer in N_samples_cancer:
                for threshold_normal_cohort in Threshold_normal_cohort:
                    for n_samples_normal in N_samples_normal:
                        if (n_samples_cancer is None) ^ (threshold_cancer_cohort is None):
                            continue
                        if (n_samples_normal is None) and (threshold_normal_cohort is None):
                            continue

                        adjusted_threshold_col = 'tmp_cancer_cohort'
                        max_threshold_col = 'tmp_normal_Nmax_sup{}'.format(threshold_normal_cohort)
                        max_threshold_col_base = 'tmp_normal_Nmax_sup{}'.format(0)


                        df = df_load.copy()
                        # Make correction for number samples passing theshold in cohort. We want to exclude the target sample in counting
                        if (n_samples_cancer is not None) and (threshold_cancer_cohort is not None):
                            df[adjusted_threshold_col] = df[get_threshold_colname(threshold_cancer_cohort, tag_cancer)]
                            df.loc[df[target_sample] >= threshold_cancer_cohort, adjusted_threshold_col] -=1 

                        # Number of kmers expressed in sample 
                        df = filter_single_col(df, 0, target_sample)
                        output_count(df, report_count, report_steps, 'Init_Sample')

                        # Number of kmers >= threshold in sample 
                        df = filter_single_col(df, threshold_target, target_sample)
                        output_count(df, report_count, report_steps, 'Filter_Sample')

                        
                        # Filter for cancer cancer cohort 
                        if (n_samples_cancer is not None) and (threshold_cancer_cohort is not None):
                            df = filter_cancer_cohort(df, n_samples_cancer, adjusted_threshold_col)
                        output_count(df, report_count, report_steps, 'Filter_Sample_Cohort')
                            


                        # Expression in gtex cohort >= threshold 
                        if threshold_normal_cohort is not None:
                            recurrence_custom =  max_recurrence_over_kmer(df, 
                                                                          get_threshold_colname(threshold_normal_cohort, tag_normal), 
                                                                          max_threshold_col)

                        # Expression in gtex cohort > 0  
                        recurrence_custom_base = max_recurrence_over_kmer(df, 
                                                                          get_threshold_colname(0.0, tag_normal),
                                                                          max_threshold_col_base) 




                        # Perform Background filtering 
                        if threshold_normal_cohort is not None:
                            df = df.merge(recurrence_custom, on = 'kmer', how = 'left')
                        df = df.merge(recurrence_custom_base, on = 'kmer', how = 'left')
                        
                        if (threshold_normal_cohort is not None) and (n_samples_normal is not None):
                            df = df.loc[ ~ ((df[max_threshold_col] >= 1) | (df[max_threshold_col_base] >= n_samples_normal)), :]
                            threshold_normal_cohort_save = threshold_normal_cohort
                            n_samples_normal_save = n_samples_normal
                        elif threshold_normal_cohort is None:
                            df = df.loc[ ~ (df[max_threshold_col_base] >= n_samples_normal), :]
                            threshold_normal_cohort_save = 'Any'
                        elif n_samples_normal is None:
                            df = df.loc[ ~ (df[max_threshold_col] >= 1), :]
                            n_samples_normal_save = 'Any'

          
                        output_count(df, report_count, report_steps, 'Filter_Sample_Cohort_CohortNormal')
                    
                        #Perform Annotated junctions filtering 
                        if filter_annot:
                            df = df[df['isAnnotated'].isna()]
                            output_count(df, report_count, report_steps, 'Filter_Sample_Cohort_CohortNormal_pepAnnot')

#                         # DEV: Exclude genes where GTEX is missing
#                         df = df.loc[df['exclude'].isna()]
#                         output_count(df, report_count, report_steps, 'Filter_Sample_Cohort_CohortNormal_pepAnnot_EXPGTEX')

                        
                        
                        
                        # Save outputs 
                        # outpaths
                        base_path_final = os.path.join(output_dir,
                                                       (f'{tag_prefix}{cancer_sample_ori}_'
                                                        f'SampleLim{threshold_target}'
                                                        f'CohortLim{threshold_cancer_cohort}'
                                                        f'Across{n_samples_cancer}_'
                                                        f'FiltNormals{save_tag}'
                                                        f'Cohortlim{threshold_normal_cohort_save}'
                                                        f'Across{n_samples_normal_save}.tsv.gz'))
                        print(f'Saving outputs to: {base_path_final} \n')
                        df.loc[:, metadata_save].to_csv(base_path_final, compression = 'gzip', index = None, sep = '\t')


    save_output_count(summary_path, report_count, report_steps, '', cancer_sample_ori, mutation_mode,
                      threshold_target, threshold_cancer_cohort, n_samples_cancer,
                          threshold_normal_cohort, n_samples_normal, save_tag)


