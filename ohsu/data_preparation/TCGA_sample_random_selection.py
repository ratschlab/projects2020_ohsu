#!/usr/bin/env python3

"""
TCGA_sample_selection.py

sample run:
time python ../junction-filtering/TCGA_sample_selection.py
-e ../immunotherapy/eth/2020_corrected_jxs/
-c ../immunotherapy/proteomics/iTRAQ_sample_mapping/

"""

import argparse
import sys
import csv; csv.field_size_limit(sys.maxsize)
from datetime import datetime
import glob
import mmh3
import os
import pandas as pd
import random
from utilities import _TCGA_ABBR_TO_CAN


def collect_cptac_samples(itraq_dir, can_main):
    filename = 'CPTAC_TCGA*{}*iTRAQ_Sample_Mapping*'.format(can_main)
    glob_form = os.path.join(itraq_dir, filename)
    itraq_file = glob.glob(glob_form)[0]
    if can_main == 'Breast':
        header_row = 0
    else:
        header_row = 11
    itraq_df = pd.read_excel(
        itraq_file, sheet_name='Proteome', header=header_row,
        usecols=['114-Biospecimen', '115-Biospecimen',
                 '116-Biospecimen', '117-Biospecimen']
    )
    itraq_df.dropna(inplace=True)
    sample_set = set()
    for col in itraq_df:
        for item in itraq_df[col]:
            try:
                if item.startswith('TCGA'):
                    sample_set.add(item[:12])
            except AttributeError:
                'encountered attribute error here:'
                print(item)
                raise AttributeError

    return sample_set


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Select random TCGA samples for junction filter analysis.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Give path to store output count/sample files.'
    )
    parser.add_argument(
        '--cptac-itraq-directory', '-c', required=True,
        help='Provide the directory containing iTRAQ sample mapping files for '
             'TCGA breast and ovarian cancer CPTAC data.'
    )
    parser.add_argument(
        '--eth-neojunctions', '-e',
        help='specify the directory containing the neojx file from '
             'doi:10.1016/j.ccell.2018.07.001: '
             'tss_complexity_counts.whitelisted.G0.01.globsum20.filtLib'
             '.conf3_neojunctions_corrected.tsv'
    )

    args = parser.parse_args()
    out_path = args.output_path
    itraq_dir = args.cptac_itraq_directory
    eth_dir = args.eth_neojunctions

    cancers=['BRCA', 'OV']
    eth_filename = (
        'tss_complexity_counts.whitelisted.G0.01.globsum20.filtLib.'
        'conf3_neojunctions_corrected.tsv'
    )
    eth_file = os.path.join(eth_dir, eth_filename)
    eth_samples = set()
    with open(eth_file) as k_input:
        k_jx_lines = csv.reader(k_input, delimiter='\t')
        for line in k_jx_lines:
            full_id, gene, jx_coords = line
            tcga_id = full_id[:12]
            eth_samples.add(tcga_id)

    for abbr in cancers:
        cancer = _TCGA_ABBR_TO_CAN[abbr]
        can_main = cancer.split('_')[0]
        cptac_samples = collect_cptac_samples(itraq_dir, can_main)
        shortlist = list(cptac_samples.intersection(eth_samples))
        shortlist.sort()
        random.seed(abs(mmh3.hash(cancer)))
        final_samples = random.sample(shortlist, 5)
        filename = '{}_5sample_random_final.csv'.format(abbr)
        with open(os.path.join(out_path, filename), 'w') as output:
            for sample in final_samples:
                output.write('{}\n'.format(sample))
                