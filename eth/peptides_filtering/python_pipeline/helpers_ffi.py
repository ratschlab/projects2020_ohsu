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
 

def filter_cancer_cohort(df, n_samples, threshold_column ):
    df = df.loc[df[threshold_column] >= n_samples]    
    return df


def get_threshold_colname(threshold, tag):
    if (threshold is not None) and (threshold > 0 ):
        threshold_column = f'{tag}filter >={threshold}'
    else:
        threshold_column = f'{tag}filter >{threshold}'
    return threshold_column


def filter_single_col(df, threshold, colname):
    if threshold:
        df = df.loc[df[colname] >= threshold, :]
    else:
        df = df.loc[df[colname] >threshold, :]
    return df


def max_recurrence_over_kmer(df, threshold_column, new_maxcol):
    df = df[['kmer', threshold_column]].groupby('kmer').max()
    df = df.reset_index()
    df = df.rename({threshold_column: new_maxcol}, axis = 1)
    return df 


def output_count(df, report_count, report_step, step_string, perform_count=True):
    '''
    Performs a count operation on the number of kmers present in spark dataframe after a given filtering step
    Note: This operation is expensive but useful if the user is interested in intermediate filtering steps
    :param perform_count: bool whether to perform a count operation
    :param df: dataframe with kmer expression counts
    :param report_count: list to store result of successive counting operations
    :param report_step: list to store name of successive counting operations
    :param step_string: str name of the counting operation
    '''
    if perform_count:
        mycount = df['kmer'].unique().shape[0]
        report_count.append(mycount)
        report_step.append(step_string)
        print(f'# {step_string} n = {mycount} kmers')


def save_output_count(output_count, report_count, report_steps, prefix, cancer_sample_ori, mutation_mode,
                      sample_expr_support_cancer, cohort_expr_support_cancer, n_samples_lim_cancer,
                          cohort_expr_support_normal, n_samples_lim_normal, id_normals):
    '''
    Saves the number of kmers present in spark dataframe after each filtering step in a tabular file
    :param output_count: str path for count file of intermediate filtering steps
    :param report_count: list to store result of successive counting operations
    :param report_step: list to store name of successive counting operations
    :param prefix: str information to be added to the result line in an info column
    :param cancer_sample_ori: str id of target cancer sample which was filtered
    :param mutation_mode: str information about whether mutations where applied or not
    :param sample_expr_support_cancer: float normalized expression threshold for the cancer target sample
    :param cohort_expr_support_cancer: float normalized expression threshold for the cancer cohort
    excluding the target sample
    hich should be met in n samples
    :param n_samples_lim_cancer: int number of cancer samples in which the cancer cohort expression threshold
    should be met
    :param cohort_expr_support_normal: float normalized expression threshold for the normal cohort
    required in any sample (>=1)
    :param n_samples_lim_normal: int number of normal samples in which any number of reads is required (>0)
    :param id_normals: str id of the normal cohort (example gtex)
    '''
    pipeline = 'peptide-based'
    if output_count:
        header = (f'{"sample"}\t{"mutation_mode"}\t{"pipeline"}\t{"min_sample_reads"}\t{"#_of_cohort_samples"}\t'
                  f'{"reads_per_cohort_sample"}\t{"#_normal_samples_allowed"}\t{"normal_cohort_id"}'
                  f'\t{"reads_per_normal_sample"}')
        line =   (f'{cancer_sample_ori}\t{mutation_mode}\t{pipeline}\t{sample_expr_support_cancer}\t{n_samples_lim_cancer}'
                  f'\t{cohort_expr_support_cancer}\t{n_samples_lim_normal}\t{id_normals}'
                  f'\t{cohort_expr_support_normal}')

        for idx in np.arange(len(report_count)):
            header += f'\t{report_steps[idx]}'
            line += f'\t{report_count[idx]}'
        if prefix:
            header += f'\t{"info"}'
            line += f'\t{prefix}'
        header += "\n"
        line += "\n"
        
#         print(header, len(header.split('\t')))
#         print(line, len(line.split('\t')))
        if not os.path.exists(output_count):
            with open(output_count,"w") as f:
                f.write(header)
        with open(output_count, "a") as f:
            f.write(line)
        logging.info(f'Save intermediate info to {output_count}')


