#!/usr/bin/env python3

"""
prep_proteomics_data.py
Python 3.6 code for preparing fasta files for CPTAC queries

SAMPLE RUN: time python ../junction-filtering/scripts/


"""
import argparse
from datetime import datetime
import glob
from matplotlib import use; use('pdf')
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import pandas as pd
import pickle
import sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from utilities import _TARGET_FILE_CANCER_MAP


def clvsites_vs_exonlength_fastas(out_dir, now):
    """

    :param jx_df:o
    :param out_dir:
    :param ref_fasta:
    :param samtools:
    :param gtf_file:
    :param uniprot_file:
    :return:
    """
    print('\nchecking cleavage sites per length from peptides')
    fasta_files = glob.glob(os.path.join(out_dir, 'fasta_info_dict*.pickle'))
    fasta_dict = {}
    for file in fasta_files:
        with open(file, 'rb') as recover:
            fasta_dict.update(pickle.load(recover))

    plot_dict = {'exon_size': [], 'cleavage_sites': []}

    npd = {}
    for junction, biexons in fasta_dict.items():
        for biexon, info in biexons.items():
            if not biexon:
                continue
            jx_pos = info['jx_pos']
            left = biexon[:jx_pos + 1]
            right = biexon[jx_pos + 1:]
            for exon_peptide in [left, right]:
                length = len(exon_peptide)
                clv_site_num = exon_peptide.count('K')
                clv_site_num += exon_peptide.count('R')
                clv_site_num -= exon_peptide.count('KP')
                clv_site_num -= exon_peptide.count('RP')
                plot_dict['exon_size'].append(length)
                plot_dict['cleavage_sites'].append(clv_site_num)
                if length in npd.keys():
                    try:
                        npd[length][clv_site_num] += 1
                    except KeyError:
                        npd[length][clv_site_num] = 1
                else:
                    npd[length] = {clv_site_num: 1}

    plot_dict = {'exon_size': [], 'cleavage_sites': [], 'count': []}
    for lengths, data in npd.items():
        for clv_site_count, count in data.items():
            plot_dict['exon_size'].append(lengths)
            plot_dict['cleavage_sites'].append(clv_site_count)
            # plot_dict['count'].append(int(count*10))
            plot_dict['count'].append(int(count))

    data_df = pd.DataFrame(plot_dict)
    print(data_df.head())
    print(data_df.shape)
    print(data_df.dtypes)

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    plt.figure(dpi=1200)
    mpl.style.use('seaborn-whitegrid')
    edgecolor = 'black'
    # color_dict = {'75to60': '#fc5a50', '90to75': '#2c6fbb', 'auc': '#020035'}
    facecolor = '#fcb001'
    plt.scatter(
        x='exon_size', y='cleavage_sites', data=data_df, c=facecolor,
        edgecolors=edgecolor,
        s='count',
        # s=10,
        linewidths=0.07, label=None,
        rasterized=True
    )
    data_df['count'] = data_df['count'].apply(lambda x: int(x / 10))
    data_df.rename({'exon_length': 'biexon_length'}, axis=1, inplace=True)
    with open(os.path.join(out_dir, 'clv_counts_peptides.tsv'), 'w') as output:
        data_df.to_csv(output, sep='\t', index=False)
    # plt.yscale('log')
    # plt.xscale('log')
    ax=plt.gca()
    ax.set_ylim(ymin=-0.5, ymax=10)
    ax.set_xlim(xmin=-5, xmax=100)
    fig_name = 'cleavage_sites_vs_exon_len_frompeptides_{}.pdf'.format(now)
    fig = plt.gcf()
    fig_file = os.path.join(out_dir, fig_name)
    fig.savefig(fig_file)
    plt.clf()
    return


def clvsites_vs_exonlength_precompute(out_dir, now):
    """

    :param jx_df:o
    :param out_dir:
    :param ref_fasta:
    :param samtools:
    :param gtf_file:
    :param uniprot_file:
    :return:
    """
    print('\nchecking cleavage sites per length from precomputed data')
    plot_files = glob.glob(os.path.join(out_dir, 'len_vs_cleavage_*.pickle'))
    plot_dict = {'exon_size': [], 'cleavage_sites': []}
    # for file in plot_files:
    #     with open(file, 'rb') as recover:
    #         temp_tryp_sites = pickle.load(recover)
    #     print(
    #         len(temp_tryp_sites['exon_size']),
    #         len(temp_tryp_sites['cleavage_sites'])
    #     )
    #     plot_dict['exon_size'].extend(temp_tryp_sites['exon_size'])
    #     plot_dict['cleavage_sites'].extend(temp_tryp_sites['cleavage_sites'])
    #
    # data_df = pd.DataFrame(plot_dict)
    npd = {}
    for file in plot_files:
        with open(file, 'rb') as recover:
            temp_tryp_sites = pickle.load(recover)
        print(
            len(temp_tryp_sites['exon_size']),
            len(temp_tryp_sites['cleavage_sites'])
        )
        plot_dict['exon_size'].extend(temp_tryp_sites['exon_size'])
        plot_dict['cleavage_sites'].extend(temp_tryp_sites['cleavage_sites'])

    len_clv = zip(plot_dict['exon_size'], plot_dict['cleavage_sites'])
    for length, clv_site_num in len_clv:
        if length in npd.keys():
            try:
                npd[length][clv_site_num] += 1
            except KeyError:
                npd[length] = {clv_site_num: 1}
        else:
            npd[length] = {clv_site_num: 1}

    # fasta_files = glob.glob(os.path.join(out_dir, 'fasta_info_dict*.pickle'))
    # fasta_dict = {}
    # for file in fasta_files:
    #     with open(file, 'rb') as recover:
    #         fasta_dict.update(pickle.load(recover))
    #
    # plot_dict = {'exon_size': [], 'cleavage_sites': []}
    #
    # npd = {}
    # for junction, biexons in fasta_dict.items():
    #     for biexon, info in biexons.items():
    #         if not biexon:
    #             continue
    #         jx_pos = info['jx_pos']
    #         left = biexon[:jx_pos + 1]
    #         right = biexon[jx_pos + 1:]
    #         for exon_peptide in [left, right]:
    #             length = len(exon_peptide)
    #             clv_site_num = exon_peptide.count('K')
    #             clv_site_num += exon_peptide.count('R')
    #             clv_site_num -= exon_peptide.count('KP')
    #             clv_site_num -= exon_peptide.count('RP')
    #             plot_dict['exon_size'].append(length)
    #             plot_dict['cleavage_sites'].append(clv_site_num)
    #             if length in npd.keys():
    #                 try:
    #                     npd[length][clv_site_num] += 1
    #                 except KeyError:
    #                     npd[length][clv_site_num] = 1
    #             else:
    #                 npd[length] = {clv_site_num: 1}

    plot_dict = {'exon_size': [], 'cleavage_sites': [], 'count': []}
    for lengths, data in npd.items():
        for clv_site_count, count in data.items():
            plot_dict['exon_size'].append(lengths)
            plot_dict['cleavage_sites'].append(clv_site_count)
            # plot_dict['count'].append(int(count*10))
            plot_dict['count'].append(int(count))

    data_df = pd.DataFrame(plot_dict)
    print(data_df.head())
    print(data_df.shape)
    print(data_df.dtypes)

    plt.rcParams.update({'figure.autolayout': True})
    plt.rcParams['figure.figsize'] = 5.0, 5.0
    plt.figure(dpi=1200)
    mpl.style.use('seaborn-whitegrid')
    edgecolor = 'black'
    # color_dict = {'75to60': '#fc5a50', '90to75': '#2c6fbb', 'auc': '#020035'}
    facecolor = '#fcb001'
    plt.scatter(
        x='exon_size', y='cleavage_sites', data=data_df, c=facecolor,
        edgecolors=edgecolor,
        s='count',
        # s=10,
        linewidths=0.07, label=None,
        rasterized=True
    )
    data_df['count'] = data_df['count'].apply(lambda x: int(x / 10))
    data_df.rename({'exon_length': 'biexon_length'}, axis=1, inplace=True)
    with open(os.path.join(out_dir, 'clv_counts_precomp.tsv'), 'w') as output:
        data_df.to_csv(output, sep='\t', index=False)
    # plt.yscale('log')
    # plt.xscale('log')
    ax = plt.gca()
    ax.set_ylim(ymin=-0.5, ymax=10)
    ax.set_xlim(xmin=-5, xmax=100)
    fig_name = 'cleavage_sites_vs_exon_len_precompute_{}.pdf'.format(now)
    fig = plt.gcf()
    fig_file = os.path.join(out_dir, fig_name)
    fig.savefig(fig_file)
    plt.clf()
    return


def create_fasta_files(fasta_dir, expts_dir, out_dir):
    """

    :param jx_df:o
    :param out_dir:
    :param ref_fasta:
    :param samtools:
    :param gtf_file:
    :param uniprot_file:
    :return:
    """

    fasta_files = glob.glob(os.path.join(fasta_dir, 'fasta_info_dict*.pickle'))
    fasta_dict = {}
    for file in fasta_files:
        with open(file, 'rb') as recover:
            fasta_dict.update(pickle.load(recover))

    fasta_line_length = 70

    for sample in _TARGET_FILE_CANCER_MAP.keys():
        expt_file = os.path.join(
            expts_dir, 'J_{}_experiments_per_peptide.tsv'.format(sample)
        )
        try:
            pep_df = pd.read_table(expt_file)
        except FileNotFoundError:
            continue
        target_peps = set(pep_df['peptide_sequence'].tolist())
        pep_ids = pep_df.set_index('peptide_sequence')['peptide_id'].to_dict()
        output_fasta = os.path.join(
            out_dir, 'J_{}_pool_kmer.fa'.format(sample)
        )
        info_to_write = set()
        for junction, biexons in fasta_dict.items():
            for pep_seq, info in biexons.items():
                if not pep_seq or pep_seq not in target_peps:
                    continue
                pep_id = pep_ids[pep_seq]
                # Add info to set: some peptides are repeated multiple times
                # but we want to ensure that each peptide-junction position
                # added is unique.
                info_to_write.add((pep_seq, (
                    pep_id, info['jx_pos'], info['between_codons'],
                    info['fiveprime_end'], info['threeprime_end'],
                    '|'.join(junction.split(';'))
                )))

        with open(output_fasta, 'w') as output:
            for pep_seq, pep_info in info_to_write:
                # peps = [
                #     pep_seq[i:i+fasta_line_length]
                #     for i in range(0, len(pep_seq), fasta_line_length)
                # ]
                peps = pep_seq
                output.write(
                    ">pepID-{};jx_pos-{};between_codons-{};"
                    "includes_5'-{};includes_3'-{};jx_coord-{}\n"
                    "".format(*pep_info)
                )
                # for pep in peps:
                #     output.write('{}\n'.format(pep))
                output.write('{}\n'.format(peps))

    return


def create_single_fasta(fasta_dir, expts_dir, out_dir):
    """

    :param jx_df:o
    :param out_dir:
    :param ref_fasta:
    :param samtools:
    :param gtf_file:
    :param uniprot_file:
    :return:
    """

    fasta_files = glob.glob(os.path.join(fasta_dir, 'fasta_info_dict*.pickle'))
    fasta_dict = {}
    for file in fasta_files:
        with open(file, 'rb') as recover:
            fasta_dict.update(pickle.load(recover))

    fasta_line_length = 70
    pep_df = pd.DataFrame()
    for sample, cancer in _TARGET_FILE_CANCER_MAP.items():
        if cancer == 'OV':
            continue
        expt_file = os.path.join(
            expts_dir, 'J_{}_experiments_per_peptide.tsv'.format(sample)
        )
        try:
            sub_df = pd.read_table(expt_file)
        except FileNotFoundError:
            continue
        pep_df = pd.concat([pep_df, sub_df])
    target_peps = set(pep_df['peptide_sequence'].tolist())
    pep_ids = pep_df.set_index('peptide_sequence')['peptide_id'].to_dict()
    output_fasta = os.path.join(out_dir, 'J_allBRCA_query_peptides.fa')
    info_to_write = set()
    for junction, biexons in fasta_dict.items():
        for pep_seq, info in biexons.items():
            if not pep_seq or pep_seq not in target_peps:
                continue
            pep_id = pep_ids[pep_seq]
            # Add info to set: some peptides are repeated multiple times
            # but we want to ensure that each peptide-junction position
            # added is unique.
            info_to_write.add((pep_seq, (
                pep_id, info['jx_pos'], info['between_codons'],
                info['fiveprime_end'], info['threeprime_end'],
                '|'.join(junction.split(';'))
            )))

    with open(output_fasta, 'w') as output:
        for pep_seq, pep_info in info_to_write:
            # peps = [
            #     pep_seq[i:i+fasta_line_length]
            #     for i in range(0, len(pep_seq), fasta_line_length)
            # ]
            peps = pep_seq
            output.write(
                ">pepID-{};jx_pos-{};between_codons-{};"
                "includes_5'-{};includes_3'-{};jx_coord-{}\n".format(*pep_info)
            )
            # for pep in peps:
            #     output.write('{}\n'.format(pep))
            output.write('{}\n'.format(peps))
    return



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collects target junctions from GTEx & TCGA SJ.out files.'
    )
    parser.add_argument(
        '--peptide-info-directory', '-f', required=True,
        help='Give the directory with the "fasta_info_dict[batch].pickle" '
             'files generated by process_shortlist_jxs.py.'
    )
    parser.add_argument(
        '--expts-per-peptide', '-e', required=True,
        help='Give the directory with "[tcga id]experiments_per_peptide.tsv" '
             'files  generated by process_shortlist_jxs.py.'
    )
    parser.add_argument(
        '--output-directory', '-o',
        help='directory in which to store output files; if no '
             'directory specified, this will be set as the same directory '
             'containing the range summarized results.'
    )

    args = parser.parse_args()
    fasta_dir = args.peptide_info_directory
    expts_dir = args.expts_per_peptide
    out_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    # clvsites_vs_exonlength_fastas(fasta_dir, now)
    # clvsites_vs_exonlength_precompute(fasta_dir, now)
    create_fasta_files(fasta_dir, expts_dir, out_dir)
    create_single_fasta(fasta_dir, expts_dir, out_dir)
    