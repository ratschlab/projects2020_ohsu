#!/usr/bin/env python3

"""
process_shortlist_jxs.py
Python 3.6 code for processing collected shortlist TCGA junctions


"""
import argparse
from datetime import datetime
import glob
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import pandas as pd
import sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from utilities import _TARGET_FILE_CANCER_MAP, _COHORT_COLORS


def plot_lines_22panel_Apr2024(j_filtered_df_dir, g_filtered_df_dir, out_dir,
                               now, logscale=False):
    j_color = '#1d8ea9'
    g_color = '#f27700'
    transparancy = 0.4
    width = 0.1
    j_positions = [1, 2, 3, 4, 5, 6]
    j_x_labels = [
        'all kmers', 'sample expression\njunction filter',
        'cancer-type cohort\njunction filter',
        'GTEx normal\njunction filter', 'splice motif\njunction filter',
        'Uniprot normal\nproteome filter'
    ]
    g_positions = [1, 2, 3, 4, 5]
    g_x_labels = [
        'Uniprot normal\nproteome filter',
        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    target_j_rows = [
        'Init_cancer', 'Filter_Sample', 'Filter_Sample_Cohort',
        'Filter_Sample_Cohort_CohortBackground', 'Filter_Motif',
        'Filter_Sample_Cohort_CohortBackground_Uniprot'
    ]
    g_cols = [11, 12, 13, 14, 15]
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.set_figheight(7)
    fig.set_figwidth(8)

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower left
            ax = axs[1, 0]
            ax.set_ylabel(
                'remaining 9-mers (#)', fontsize=12
            )
            ax.set_xticks(j_positions,)
            ax.set_xticklabels(
                j_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )

        else:
            # upper left
            ax = axs[0, 0]
            ax.set_xticks(j_positions)
            ax.set_xticklabels(
                ['', '', '', '', '', ''],
            )
            ax.set_title('JP', fontsize=18)
            ax.set_ylabel(
                'remaining 9-mers (#)', fontsize=12
            )
        color = j_color
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.EngFormatter('', sep=''))
        j_file = glob.glob(os.path.join(
            j_filtered_df_dir, 'J_filtered_df_{}*.tsv'.format(samp)
        ))
        j_samp_df = pd.read_table(j_file[0], sep='\t')
        for _, row in j_samp_df[target_j_rows].iterrows():
            ax.plot(
                j_positions, row, color=color, alpha=transparancy,
                linewidth=width
            )

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower right
            ax = axs[1, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('OV', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_xticklabels(
                g_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )
        else:
            # upper right
            ax = axs[0, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('BRCA', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_title('GP', fontsize=18)
            ax.set_xticklabels(['', '', '', '', ''])
        color = g_color
        ax.tick_params(axis='both', which='major', labelsize=10)

        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.EngFormatter('', sep=''))
        g_file = glob.glob(os.path.join(
            g_filtered_df_dir, 'G_filtered_df_{}*.tsv'.format(samp)
        ))
        g_samp_df = pd.read_table(
            g_file[0], usecols=g_cols, names=g_x_labels, sep='\t'
        )
        for _, row in g_samp_df.iterrows():
            ax.plot(
                g_positions, row, color=color, alpha=transparancy,
                linewidth=width
            )

    ax = axs[0, 0]
    if logscale:
        ax.set_yscale('log')
        fig_name = 'lineplots_remaining_count_log_{}.pdf'.format(now)
    else:
        fig_name = 'lineplots_remaining_count_{}.pdf'.format(now)
    fig = plt.gcf()
    fig_file = os.path.join(out_dir, fig_name)
    fig.savefig(
        fig_file,
        bbox_inches='tight',
        pad_inches=0.1
    )
    return



def plot_fulldiff_22panel_Apr2024(j_filtered_df_dir, g_filtered_df_dir,
                                  out_dir, now, logscale=False):
    j_color = '#1d8ea9'
    g_color = '#f27700'
    transparancy = 0.4
    width = 0.1
    j_positions = [1, 2, 3, 4, 5]

    j_x_labels = [
        'sample expression\njunction filter',
        'cancer-type cohort\njunction filter',
        'GTEx normal\njunction filter', 'splice motif\njunction filter',
        'Uniprot normal\nproteome filter'
    ]
    j_total_col = 'Init_cancer'
    g_positions = [1, 2, 3, 4]
    g_x_labels = [

        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    target_j_rows = [
        'Init_cancer', 'Filter_Sample', 'Filter_Sample_Cohort',
        'Filter_Sample_Cohort_CohortBackground', 'Filter_Motif',
        'Filter_Sample_Cohort_CohortBackground_Uniprot'
    ]
    target_g_rows = [
        'Uniprot normal\nproteome filter',
        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    g_total_col = 'Uniprot normal\nproteome filter'

    g_cols = [11, 12, 13, 14, 15]
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.set_figheight(7)
    fig.set_figwidth(8)

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower left
            ax = axs[1, 0]
            ax.set_ylabel(
                'remaining 9-mers (%)', fontsize=12
            )
            ax.set_xticks(j_positions,)
            ax.set_xticklabels(
                j_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )

        else:
            # upper left
            ax = axs[0, 0]
            ax.set_xticks(j_positions)
            ax.set_xticklabels(
                ['', '', '', '', '', ''],
            )
            ax.set_title('JP', fontsize=18)
            ax.set_ylabel(
                'remaining 9-mers (%)', fontsize=12
            )
        color = j_color
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.PercentFormatter())
        j_file = glob.glob(os.path.join(
            j_filtered_df_dir, 'J_filtered_df_{}*.tsv'.format(samp)
        ))
        j_samp_df = pd.read_table(j_file[0], sep='\t')
        for _, row in j_samp_df[target_j_rows].iterrows():
            new_vals = []
            init_val = row[j_total_col]
            for j, value in enumerate(row):
                if j == 0:
                    continue
                new_vals.append(100 * value / init_val)
            ax.plot(
                j_positions, new_vals, color=color, alpha=transparancy,
                linewidth=width
            )

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower right
            ax = axs[1, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('OV', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_xticklabels(
                g_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )
        else:
            # upper right
            ax = axs[0, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('BRCA', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_title('GP', fontsize=18)
            ax.set_xticklabels(['', '', '', '', ''])
        color = g_color
        ax.tick_params(axis='both', which='major', labelsize=10)

        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.PercentFormatter())
        g_file = glob.glob(os.path.join(
            g_filtered_df_dir, 'G_filtered_df_{}*.tsv'.format(samp)
        ))
        g_samp_df = pd.read_table(
            g_file[0], usecols=g_cols, names=target_g_rows, sep='\t'
        )
        for _, row in g_samp_df.iterrows():
            new_vals = []
            init_val = row[g_total_col]
            for j, value in enumerate(row):
                if j == 0:
                    continue
                new_vals.append(100 * value / init_val)
            ax.plot(
                g_positions, new_vals, color=color, alpha=transparancy,
                linewidth=width
            )

    ax = axs[0, 0]
    if logscale:
        ax.set_yscale('log')
        fig_name = 'lineplots_difference_from_init_log_{}.pdf'.format(now)
    else:
        fig_name = 'lineplots_difference_from_init_{}.pdf'.format(now)
    fig = plt.gcf()
    fig_file = os.path.join(out_dir, fig_name)
    fig.savefig(
        fig_file,
        bbox_inches='tight',
        pad_inches=0.1
    )
    return


def plot_stepdiff_22panel_Apr2024(j_filtered_df_dir, g_filtered_df_dir,
                                  out_dir, now, logscale=False):
    j_color = '#1d8ea9'
    g_color = '#f27700'
    transparancy = 0.4
    width = 0.1
    j_positions = [1, 2, 3, 4, 5]
    j_x_labels = [
        'sample expression\njunction filter',
        'cancer-type cohort\njunction filter',
        'GTEx normal\njunction filter', 'splice motif\njunction filter',
        'Uniprot normal\nproteome filter'
    ]
    j_total_col = 'Init_cancer'
    g_positions = [1, 2, 3, 4]
    g_x_labels = [

        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    target_j_rows = [
        'Init_cancer', 'Filter_Sample', 'Filter_Sample_Cohort',
        'Filter_Sample_Cohort_CohortBackground', 'Filter_Motif',
        'Filter_Sample_Cohort_CohortBackground_Uniprot'
    ]
    target_g_rows = [
        'Uniprot normal\nproteome filter',
        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    g_total_col = 'Uniprot normal\nproteome filter'

    g_cols = [11, 12, 13, 14, 15]
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.set_figheight(7)
    fig.set_figwidth(8)

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower left
            ax = axs[1, 0]
            ax.set_ylabel(
                'remaining 9-mers\n(% of previous step)', fontsize=12
            )
            ax.set_xticks(j_positions,)
            ax.set_xticklabels(
                j_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )

        else:
            # upper left
            ax = axs[0, 0]
            ax.set_xticks(j_positions)
            ax.set_xticklabels(
                ['', '', '', '', '', ''],
            )
            ax.set_title('JP', fontsize=18)
            ax.set_ylabel(
                'remaining 9-mers\n(% of previous step)', fontsize=12
            )
        color = j_color
        ax.tick_params(axis='both', which='major', labelsize=10)
        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.PercentFormatter())
        j_file = glob.glob(os.path.join(
            j_filtered_df_dir, 'J_filtered_df_{}*.tsv'.format(samp)
        ))
        j_samp_df = pd.read_table(j_file[0], sep='\t')
        for _, row in j_samp_df[target_j_rows].iterrows():
            new_vals = []
            prev_val = row[j_total_col]
            for j, value in enumerate(row):
                if j == 0:
                    continue
                try:
                    new_vals.append(100 * value / prev_val)
                except ZeroDivisionError:
                    new_vals.append(100)
                prev_val = value
            ax.plot(
                j_positions, new_vals, color=color, alpha=transparancy,
                linewidth=width
            )

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower right
            ax = axs[1, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('OV', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_xticklabels(
                g_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )
        else:
            # upper right
            ax = axs[0, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('BRCA', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_title('GP', fontsize=18)
            ax.set_xticklabels(['', '', '', '', ''])
        color = g_color
        ax.tick_params(axis='both', which='major', labelsize=10)

        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.PercentFormatter())
        g_file = glob.glob(os.path.join(
            g_filtered_df_dir, 'G_filtered_df_{}*.tsv'.format(samp)
        ))
        g_samp_df = pd.read_table(
            g_file[0], usecols=g_cols, names=target_g_rows, sep='\t'
        )
        for _, row in g_samp_df.iterrows():
            new_vals = []
            prev_val = row[g_total_col]
            for j, value in enumerate(row):
                if j == 0:
                    continue
                try:
                    new_vals.append(100 * value / prev_val)
                except ZeroDivisionError:
                    new_vals.append(100)
                prev_val = value
            ax.plot(
                g_positions, new_vals, color=color, alpha=transparancy,
                linewidth=width
            )

    ax = axs[0, 0]
    if logscale:
        ax.set_yscale('log')
        fig_name = 'lineplots_difference_from_prev_log_{}.pdf'.format(now)
    else:
        fig_name = 'lineplots_difference_from_prev_{}.pdf'.format(now)
    fig = plt.gcf()
    fig_file = os.path.join(out_dir, fig_name)
    fig.savefig(
        fig_file,
        bbox_inches='tight',
        pad_inches=0.1
    )
    return


def plot_stepdiff_firstincl_22panel_Apr2024(j_filtered_df_dir,
                                            g_filtered_df_dir, out_dir, now,
                                            logscale=False):
    j_color = '#1d8ea9'
    g_color = '#f27700'
    transparancy = 0.4
    width = 0.1
    j_positions = [1, 2, 3, 4, 5, 6]
    j_x_labels = [
        'initial sample',
        'sample expression\njunction filter',
        'cancer-type cohort\njunction filter',
        'GTEx normal\njunction filter', 'splice motif\njunction filter',
        'Uniprot normal\nproteome filter'
    ]
    j_total_col = 'Init_cancer'
    g_positions = [1, 2, 3, 4, 5]
    g_x_labels = [
        'Uniprot normal\nproteome filter',
        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    target_j_rows = [
        'Init_cancer', 'Filter_Sample', 'Filter_Sample_Cohort',
        'Filter_Sample_Cohort_CohortBackground', 'Filter_Motif',
        'Filter_Sample_Cohort_CohortBackground_Uniprot'
    ]
    target_g_rows = [
        'Uniprot normal\nproteome filter',
        'sample expression\nproteome filter',
        'cancer-type cohort\nproteome filter',
        'GTEx normal\nproteome filter', 'GTEx normal\njunction filter'
    ]
    g_total_col = 'Uniprot normal\nproteome filter'

    g_cols = [11, 12, 13, 14, 15]
    fig, axs = plt.subplots(2, 2)
    plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.set_figheight(7)
    fig.set_figwidth(8)

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower left
            ax = axs[1, 0]
            ax.set_ylabel(
                'remaining 9-mers\n(% of previous step)', fontsize=12
            )
            ax.set_xticks(j_positions,)
            ax.set_xticklabels(
                j_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )

        else:
            # upper left
            ax = axs[0, 0]
            ax.set_xticks(j_positions)
            ax.set_xticklabels(
                ['', '', '', '', '', ''],
            )
            ax.set_title('JP', fontsize=18)
            ax.set_ylabel(
                'remaining 9-mers\n(% of previous step)', fontsize=12
            )
        color = j_color
        ax.tick_params(axis='both', which='major', labelsize=10)

        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.PercentFormatter())
        j_file = glob.glob(os.path.join(
            j_filtered_df_dir, 'J_filtered_df_{}*.tsv'.format(samp)
        ))
        j_samp_df = pd.read_table(j_file[0], sep='\t')
        expt_ct = 0
        for _, row in j_samp_df[target_j_rows].iterrows():
            expt_ct += 1
            new_vals = [100]
            prev_val = row[j_total_col]
            for j, value in enumerate(row):
                if j == 0:
                    continue
                try:
                    new_vals.append(100 * value / prev_val)
                except ZeroDivisionError:
                    new_vals.append(100)
                prev_val = value
            ax.plot(
                j_positions, new_vals, color=color, alpha=transparancy,
                linewidth=width
            )
        print('junctions: {} experiments, sample {}'.format(expt_ct, samp))

    for samp, can in _TARGET_FILE_CANCER_MAP.items():
        if can == 'OV':
            # lower right
            ax = axs[1, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('OV', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_xticklabels(
                g_x_labels,
                rotation=45,
                ha="right",
                fontsize=12
            )
        else:
            # upper right
            ax = axs[0, 1]
            ax.set_xticks(g_positions)
            ax2 = ax.twinx()
            ax2.set_ylabel('BRCA', fontsize=18)
            ax2.set_yticklabels([])
            ax.set_title('GP', fontsize=18)
            ax.set_xticklabels(['', '', '', '', ''])
        color = g_color
        ax.tick_params(axis='both', which='major', labelsize=10)

        ax.yaxis.grid(True)
        ax.xaxis.grid(True)
        ax.set_axisbelow(True)
        if logscale:
            ax.set_yscale('log')
        else:
            ax.yaxis.set_major_formatter(ticker.PercentFormatter())
        g_file = glob.glob(os.path.join(
            g_filtered_df_dir, 'G_filtered_df_{}*.tsv'.format(samp)
        ))
        g_samp_df = pd.read_table(
            g_file[0], usecols=g_cols, names=target_g_rows, sep='\t'
        )
        expt_ct = 0
        for _, row in g_samp_df.iterrows():
            new_vals = [100]
            prev_val = row[g_total_col]
            for j, value in enumerate(row):
                if j == 0:
                    continue
                try:
                    new_vals.append(100 * value / prev_val)
                except ZeroDivisionError:
                    new_vals.append(100)
                prev_val = value
            ax.plot(
                g_positions, new_vals, color=color, alpha=transparancy,
                linewidth=width
            )
            expt_ct += 1
        print('graph: {} experiments, sample {}'.format(expt_ct, samp))

    ax = axs[0, 0]

    if logscale:
        ax.set_yscale('log')
        fig_name = (
            'lineplots_diff_from_prev_startincluded_log_{}.pdf'.format(now)
        )
    else:
        fig_name = 'lineplots_diff_from_prev_startincluded_{}.pdf'.format(now)
    fig = plt.gcf()
    fig_file = os.path.join(out_dir, fig_name)
    fig.savefig(
        fig_file,
        bbox_inches='tight',
        pad_inches=0.1
    )
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collects target junctions from GTEx & TCGA SJ.out files.'
    )
    parser.add_argument(
        '--j-filtered-output-dir', '-j',
        help='Give the directory with filtered_df files generated by '
             'filter_shortlist_jxs.py, containing junction and peptide counts '
             'from target samples.'
    )
    parser.add_argument(
        '--g-filtered-output-dir', '-g',
        help='Give the directory with filtered_df files from the graph '
             'pipeline, containing junction and peptide counts from target '
             'samples.'
    )
    parser.add_argument(
        '--g-filtered-output-old-dir', '-G',
        help='Give the directory with filtered_df files from the graph '
             'pipeline from Nov.2021, containing junction and peptide counts '
             'from target samples.'
    )
    parser.add_argument(
        '--subset-df-dir', '-s',
        help='Give the directory with the "sub_df_[TCGA-ID].tsv" files '
             'generated by process_shortlist_jxs.py, containing collected '
             'junctions from target samples w/ aggregated normal sample '
             'counts at varying coverage cutoffs.'
    )
    parser.add_argument(
        '--output-directory', '-o', default='.',
        help='Specify separate output directory, if desired.'
    )

    args = parser.parse_args()
    j_filtered_df_dir = args.j_filtered_output_dir
    g_filtered_df_dir = args.g_filtered_output_dir
    g_old_dir = args.g_filtered_output_old_dir
    subset_dir = args.subset_df_dir
    out_dir = args.output_directory

    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')

    plot_lines_22panel_Apr2024(
        j_filtered_df_dir, g_filtered_df_dir, out_dir, now
    )
    plot_lines_22panel_Apr2024(
        j_filtered_df_dir, g_filtered_df_dir, out_dir, now, logscale=True
    )
    plot_fulldiff_22panel_Apr2024(
        j_filtered_df_dir, g_filtered_df_dir, out_dir, now, logscale=False
    )
    plot_fulldiff_22panel_Apr2024(
        j_filtered_df_dir, g_filtered_df_dir, out_dir, now, logscale=True
    )
    plot_stepdiff_22panel_Apr2024(
        j_filtered_df_dir, g_filtered_df_dir, out_dir, now
    )
    plot_stepdiff_firstincl_22panel_Apr2024(
        j_filtered_df_dir, g_filtered_df_dir, out_dir, now
    )
