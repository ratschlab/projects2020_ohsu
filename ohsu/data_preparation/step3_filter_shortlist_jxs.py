#!/usr/bin/env python3

"""
filter_shortlist_jxs.py
Python 3.6 code for filtering collected shortlist TCGA junctions


"""
import argparse
from datetime import datetime
import os
import pandas as pd
import sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from utilities import _TARGET_FILE_CANCER_MAP, _ANNOTATED_S_FILE, _DF_FILE_TAGS
from utilities import _TCGA_ALL, _GTEX_CORE, _GTEX_BR, _GTEX_OV, _GTEX_RESTRICT
from utilities import over_x_count_key, normal_filter, _MOTIF_FILTER
from utilities import _MAX_NORMAL_SAMPS, _MAX_NORMAL_READS, _NORM_FILTER
from utilities import _NF_OPTS, _MNS_OPTS, _MOTIF_OPTS, _MSR_OPTS
from utilities import _MIN_CHRT_SAMP_READS, _MIN_COHORT_SAMPS, _MNR_OPTS
from utilities import _CANON_MOTIFS, _MCS_OPTS, _MCSR_OPTS, _NORMS, _COHORT
from utilities import _ASN_COUNT_COLS, _PRE_IF_EP_COUNT, _MIN_SAMP_READS
from utilities import _SHORTLIST_FILE, _IF_BIEX, _UNIPROT_FINAL_COL
from utilities import _FILTERED_IF_EPS, _IF_PRE_EPS, _GTEX_ALL
from utilities import _FILTERED_IF_COUNT


_MOTIF = 'motif'
_STAT_FILE = 'filtered_summary_stats.tsv'
_UNDER = 'under_column'
_OVER = 'over_column'
_JX_COUNT = 'final_jx_count'
_FA_COUNT = 'frame-agnostic_pep_count'
_IF_COUNT = 'in-frame_pep_count'
_NORM_FLAG_MAPPER = {
    # 'core_GTEx': 'Gtexcore',
    'GTEx': 'GtexCohort',
    # 'core_GTEx': 'GtexCohort',
    'paired': 'Matched',
    'All': 'GtexTcga'
}


def create_filter_set_dict(cancertype, outfile):
    """Assembles all filter sets for junction filtering.

    Defines:
    - set of normal samples to use
    - whether or not to filter on motif
    - minimum number of in-sample reads required
    - minimum number of cohort reads required
    - maximum number of samples in normal cohort allowed
    - maximum number of reads per normal sample allowed

    Returns dict of assembled filters.
    """
    filter_sets = {}
    expt_info = {
        'id': [], _MOTIF_FILTER: [], _MIN_SAMP_READS: [],
        _MIN_COHORT_SAMPS: [], _MIN_CHRT_SAMP_READS: [],
        _NORM_FILTER: [], _MAX_NORMAL_SAMPS: [], _MAX_NORMAL_READS: []
    }
    expt_counter_map = {
        0: 'N', -1: 'A', 10: 'X',
        'paired': 'P', 'GTEx': 'G', 'core_GTEx': 'C', 'All': 'A',
    }
    fground_expt_counter_map = {0: 'A'}

    motif_map = {
        1: 'C', 0: 'A'
    }
    for nf in _NF_OPTS:
        for mf in _MOTIF_OPTS:
            for msr in _MSR_OPTS:
                cohort_samps = _MCS_OPTS[cancertype]
                for mcs in cohort_samps:
                    for mcr in _MCSR_OPTS:
                        if mcs == 0 and mcr > 0:
                            continue
                        # Add single filter set w/ no normal reads allowed:
                        # 0 normal samples with under 1 read.
                        # expt_counter += 1
                        expt_counter = (
                            f'J{msr}'
                            f'{fground_expt_counter_map.get(mcr, mcr)}'
                            f'{expt_counter_map.get(mcs, mcs)}'
                            f'01'
                            f'{expt_counter_map[nf]}'
                            f'{motif_map[mf]}'
                        )
                        expt_id = str(expt_counter)
                        filter_sets[expt_id] = {}
                        filter_sets[expt_id][_MOTIF_FILTER] = mf
                        filter_sets[expt_id][_MIN_SAMP_READS] = msr
                        filter_sets[expt_id][_MIN_CHRT_SAMP_READS] = mcr
                        filter_sets[expt_id][_MAX_NORMAL_READS] = 0
                        filter_sets[expt_id][_MAX_NORMAL_SAMPS] = 1
                        filter_sets[expt_id][_NORM_FILTER] = nf
                        filter_sets[expt_id][_MIN_COHORT_SAMPS] = mcs
                        expt_info['id'].append(expt_id)
                        expt_info[_MOTIF_FILTER].append(mf)
                        expt_info[_MIN_SAMP_READS].append(msr)
                        expt_info[_MIN_COHORT_SAMPS].append(mcs)
                        expt_info[_MIN_CHRT_SAMP_READS].append(mcr)
                        expt_info[_NORM_FILTER].append(nf)
                        expt_info[_MAX_NORMAL_SAMPS].append(1)
                        expt_info[_MAX_NORMAL_READS].append(0)
                        for mns in _MNS_OPTS + [-1]:
                            for mnr in _MNR_OPTS + [-1]:
                                if mns == -1 and mnr == -1:
                                    continue
                                if mns == 1:
                                    continue
                                if mnr == 0:
                                    continue
                                expt_counter = (
                                    f'J{msr}'
                                    f'{fground_expt_counter_map.get(mcr, mcr)}'
                                    f'{expt_counter_map.get(mcs, mcs)}'
                                    f'{expt_counter_map.get(mnr, mnr)}'
                                    f'{expt_counter_map.get(mns, mns)}'
                                    f'{expt_counter_map[nf]}'
                                    f'{motif_map[mf]}'
                                )
                                expt_id = str(expt_counter)
                                filter_sets[expt_id] = {}
                                filter_sets[expt_id][_MOTIF_FILTER] = mf
                                filter_sets[expt_id][
                                    _MIN_SAMP_READS] = msr
                                filter_sets[expt_id][
                                    _MIN_CHRT_SAMP_READS] = mcr
                                filter_sets[expt_id][_MAX_NORMAL_READS] = mnr
                                filter_sets[expt_id][_MAX_NORMAL_SAMPS] = mns
                                filter_sets[expt_id][_NORM_FILTER] = nf
                                filter_sets[expt_id][_MIN_COHORT_SAMPS] = mcs
                                expt_info['id'].append(expt_id)
                                expt_info[_MOTIF_FILTER].append(mf)
                                expt_info[_MIN_SAMP_READS].append(msr)
                                expt_info[_MIN_COHORT_SAMPS].append(mcs)
                                expt_info[_MIN_CHRT_SAMP_READS].append(mcr)
                                expt_info[_NORM_FILTER].append(nf)
                                expt_info[_MAX_NORMAL_SAMPS].append(mns)
                                expt_info[_MAX_NORMAL_READS].append(mnr)
    outdf = pd.DataFrame(expt_info)
    with open(outfile, 'w') as output:
        outdf.to_csv(output, index=False, sep='\t')
    return filter_sets


def apply_cohort_filters(filtered_df, expt, cancer):
    samp_count = expt[_MIN_COHORT_SAMPS] + 1
    read_count = expt[_MIN_CHRT_SAMP_READS]
    col = '{}_over_{}'.format(cancer, read_count)
    filtered_df = filtered_df.loc[filtered_df[col] >= samp_count]
    return filtered_df


def apply_normal_filters(filtered_df, expt, cancer):
    norm_filt = expt[_NORM_FILTER]
    max_norm_smps = expt[_MAX_NORMAL_SAMPS]
    max_norm_rds = expt[_MAX_NORMAL_READS]
    under_val = max(max_norm_rds, min(_MNR_OPTS))
    temp_col = 'temp_col'
    if norm_filt == 'All':
        col = [_GTEX_CORE, _GTEX_RESTRICT, _TCGA_ALL]
        if max_norm_smps == -1:
            def normal_filter_builder(row):
                return normal_filter(
                    0,
                    row[over_x_count_key(col[0], under_val)],
                    # Setting this arbitrarily: we want only to include the
                    # reads filter here, not the number of samples.
                    1,
                    0,
                    row[over_x_count_key(col[1], under_val)],
                    0,
                    row[over_x_count_key(col[2], under_val)]
                )
        elif max_norm_rds == -1:
            def normal_filter_builder(row):
                return normal_filter(
                    row[over_x_count_key(col[0], 0)],
                    0,
                    max_norm_smps,
                    row[over_x_count_key(col[1], 0)],
                    0,
                    row[over_x_count_key(col[2], 0)],
                    0
                )
        else:
            def normal_filter_builder(row):
                return normal_filter(
                    row[over_x_count_key(col[0], 0)],
                    row[over_x_count_key(col[0], under_val)],
                    max_norm_smps,
                    row[over_x_count_key(col[1], 0)],
                    row[over_x_count_key(col[1], under_val)],
                    row[over_x_count_key(col[2], 0)],
                    row[over_x_count_key(col[2], under_val)]
                )
    elif norm_filt == 'GTEx':
        col = [_GTEX_CORE, _GTEX_RESTRICT]
        if max_norm_smps == -1:
            def normal_filter_builder(row):
                return normal_filter(
                    0,
                    row[over_x_count_key(col[0], under_val)],
                    # Setting this arbitrarily: we want only to include the
                    # reads filter here, not the number of samples.
                    1,
                    0,
                    row[over_x_count_key(col[1], under_val)]
                )
        elif max_norm_rds == -1:
            def normal_filter_builder(row):
                return normal_filter(
                    row[over_x_count_key(col[0], 0)],
                    0,
                    max_norm_smps,
                    row[over_x_count_key(col[1], 0)],
                    0
                )
        else:
            def normal_filter_builder(row):
                return normal_filter(
                    row[over_x_count_key(col[0], 0)],
                    row[over_x_count_key(col[0], under_val)],
                    max_norm_smps,
                    row[over_x_count_key(col[1], 0)],
                    row[over_x_count_key(col[1], under_val)]
                )
    else:
        if norm_filt == 'core_GTEx':
            col = _GTEX_CORE
        elif cancer == 'OV':
            col = _GTEX_OV
        else:
            col = _GTEX_BR
        if max_norm_smps == -1:
            def normal_filter_builder(row):
                return normal_filter(
                    0,
                    row[over_x_count_key(col, under_val)],
                    # Setting this arbitrarily: we want only to include the
                    # reads filter here, not the number of samples.
                    1
                )
        elif max_norm_rds == -1:
            def normal_filter_builder(row):
                return normal_filter(
                    row[over_x_count_key(col, 0)],
                    0,
                    max_norm_smps
                )
        else:
            def normal_filter_builder(row):
                return normal_filter(
                    row[over_x_count_key(col, 0)],
                    row[over_x_count_key(col, under_val)],
                    max_norm_smps
                )
    filtered_df[temp_col] = filtered_df.apply(normal_filter_builder, axis=1)
    filtered_df = filtered_df.loc[filtered_df[temp_col] == 1]
    filtered_df.drop([temp_col], axis=1, inplace=True)
    return filtered_df


def filter_jxs(jx_df, out_dir, batch=None, genelist=''):
    peptide_dict = {}
    pep_ids = {}
    curr_pep_id = 100000
    jx_df.dropna(subset=[_IF_PRE_EPS], axis=0, inplace=True)
    jx_df[_FILTERED_IF_EPS].fillna('', inplace=True)
    pref_split = 'split_prefiltered'
    post_split = 'split_postfiltered'
    jx_df[pref_split] = jx_df[_IF_PRE_EPS].apply(lambda x: x.split(';'))
    pre_pep = jx_df.explode(pref_split)[['jx', pref_split]]
    pre_pep[_IF_PRE_EPS] = pre_pep[pref_split]
    jx_df[post_split] = jx_df[_FILTERED_IF_EPS].apply(lambda x: x.split(';'))
    filtpep = jx_df.explode(post_split)[['jx', post_split]]
    filtpep[_FILTERED_IF_EPS] = filtpep[post_split]
    if genelist:
        genedf = pd.read_table(genelist, header=None, names=['genes'])
        gene_set = set(genedf['genes'].tolist())
        jx_df['gene_id'].fillna('', inplace=True)
        jx_df['target_gene'] = jx_df['gene_id'].apply(
            lambda x: int(len(set(x.split(',')).intersection(gene_set)) > 0)
        )
        jx_df = jx_df.loc[jx_df['target_gene'] == 1].copy()
    plot_order = {
        _MIN_SAMP_READS: [],
        _COHORT: [_MIN_COHORT_SAMPS, _MIN_CHRT_SAMP_READS],
        _NORMS: [_NORM_FILTER, _MAX_NORMAL_SAMPS, _MAX_NORMAL_READS],
        _MOTIF_FILTER: [],
    }
    plot_cols = [
        'sample', 'mutation_mode', 'pipeline', 'expt_id', _MIN_SAMP_READS,
        _MIN_COHORT_SAMPS, _MIN_CHRT_SAMP_READS, _MAX_NORMAL_SAMPS,
        _NORM_FILTER, _MAX_NORMAL_READS, _MOTIF_FILTER, 'Init_cancer',
        _ASN_COUNT_COLS[_MIN_SAMP_READS], _ASN_COUNT_COLS[_COHORT],
        _ASN_COUNT_COLS[_NORMS], _ASN_COUNT_COLS[_MOTIF_FILTER],
        _UNIPROT_FINAL_COL
    ]
    samps_to_check = [
        'TCGA-AO-A0JM-01A-21R-A056-07',
        'TCGA-BH-A18V-01A-11R-A12D-07',
        'TCGA-BH-A18V-06A-11R-A213-07',
        'TCGA-A2-A0D2-01A-21R-A034-07',
        'TCGA-A2-A0SX-01A-12R-A084-07',
        'TCGA-C8-A12P-01A-11R-A115-07',
        'TCGA-25-1319-01A-01R-1565-13',
        'TCGA-25-1313-01A-01R-1565-13',
        'TCGA-61-2008-01A-02R-1568-13',
        'TCGA-24-1431-01A-01R-1566-13',
        'TCGA-24-2298-01A-01R-1569-13',
    ]
    if batch is not None:
        samps_to_check = samps_to_check[batch]
    for samp, cancer in _TARGET_FILE_CANCER_MAP.items():
        max_jxs = set()
        if samp not in samps_to_check:
            continue
        filename = 'J_filtered_df_{}'.format(samp)
        for filter in plot_order:
            filename += '_' + _DF_FILE_TAGS[filter]
        filename += '.tsv'
        filepath = os.path.join(out_dir, filename)
        if os.path.isfile(filepath):
            continue
        try:
            sample_df = jx_df.loc[jx_df[samp] > 0].copy()
        except KeyError:
            continue
        print('sample {}'.format(samp))
        print('{} junctions'.format(len(sample_df)))
        sample_df = sample_df.loc[sample_df[_PRE_IF_EP_COUNT] > 0].copy()
        tot_asns = len(set(
            pre_pep.loc[pre_pep['jx'].isin(sample_df['jx'])][_IF_PRE_EPS]
        ))
        plot_dict = {col: [] for col in plot_cols}
        exptfilterfile = os.path.join(
            out_dir, 'J_{}_experiment_map.tsv'.format(samp)
        )
        expt_dict = create_filter_set_dict(cancer, exptfilterfile)
        for ex_id, expt in expt_dict.items():
            plot_dict['sample'] = samp[:19]
            plot_dict['mutation_mode'] = 'ref'
            plot_dict['pipeline'] = 'junction-based'
            plot_dict['expt_id'].append(ex_id)
            sub_df = sample_df.copy()
            plot_dict['Init_cancer'].append(tot_asns)
            for overall_filter, subfilters in plot_order.items():
                if overall_filter in {_COHORT, _NORMS}:
                    for subfilt in subfilters:
                        plot_dict[subfilt].append(expt[subfilt])
                    if overall_filter == _COHORT:
                        sub_df = apply_cohort_filters(sub_df, expt, cancer)
                    elif overall_filter == _NORMS:
                        sub_df = apply_normal_filters(sub_df, expt, cancer)
                else:
                    param = expt[overall_filter]
                    plot_dict[overall_filter].append(param)
                    if overall_filter == _MOTIF_FILTER:
                        if param:
                            sub_df = sub_df.loc[
                                sub_df[_MOTIF].isin(_CANON_MOTIFS)
                            ]
                    else:
                        # overall_filter == _MIN_SAMP_READS
                        sub_df = sub_df.loc[sub_df[samp] >= param]
                plot_dict[_ASN_COUNT_COLS[overall_filter]].append(len(set(
                    pre_pep.loc[pre_pep['jx'].isin(sub_df['jx'])][_IF_PRE_EPS]
                )))
            plot_dict[_UNIPROT_FINAL_COL].append(len(set(
                filtpep.loc[filtpep['jx'].isin(sub_df['jx'])][_FILTERED_IF_EPS]
            )))
            kmer_file = f'J_{samp}_{ex_id}.tsv'
            max_jxs.update(sub_df['jx'].tolist())
            sub_df = sub_df.loc[sub_df[_FILTERED_IF_COUNT] > 0]
            with open(os.path.join(out_dir, kmer_file), 'w') as output:
                output.write('kmer\n')
                for pepstr, jx in zip(sub_df[_FILTERED_IF_EPS], sub_df['jx']):
                    try:
                        peplist = pepstr.split(';')
                    except AttributeError:
                        continue
                    for pep in peplist:
                        if pep == '':
                            continue
                        output.write(f'{pep}\t{jx}\n')

            for pepstr in sub_df[_IF_BIEX]:
                try:
                    peplist = pepstr.split(';')
                except AttributeError:
                    continue
                for pep in peplist:
                    if pep == '':
                        continue
                    if pep in pep_ids:
                        peptide_dict[pep_ids[pep]]['experiments'].append(ex_id)
                    else:
                        curr_pep_id += 1
                        pep_ids[pep] = curr_pep_id
                        peptide_dict[curr_pep_id] = {
                            'seq': pep, 'experiments': [ex_id]
                        }
        subdf_file = os.path.join(out_dir, 'sub_df_{}.tsv'.format(samp))
        sample_df = sample_df.loc[sample_df['jx'].isin(max_jxs)]
        sample_df.to_csv(subdf_file, index=False, sep='\t')
        plot_df = pd.DataFrame(plot_dict)
        plot_df.to_csv(filepath, index=False, sep='\t')
        peptide_expt_outfile = os.path.join(
            out_dir, 'J_{}_experiments_per_peptide.tsv'.format(samp)
        )
        with open(peptide_expt_outfile, 'w') as output:
            output.write('peptide_id\tpeptide_sequence\texperiment_ids\n')
            for pep_id, info in peptide_dict.items():
                output.write('{}\t{}\t{}\n'.format(
                    pep_id, info['seq'], ';'.join(info['experiments'])
                ))
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collects target junctions from GTEx & TCGA SJ.out files.'
    )
    parser.add_argument(
        '--target-junction-directory', '-j', required=True,
        help='Give the directory with the "annotated_all_sample_shortlist.tsv"'
             ' file generated by process_shortlist_jxs.py, containing '
             'collected junctions from target samples w/ aggregated normal '
             'sample counts at varying coverage cutoffs.'
    )
    parser.add_argument(
        '--batch-number', type=int,
        help='For running in batches on slurm: the slurm array task id to '
             'process correct chunk of the dataframe.'
    )
    parser.add_argument(
        '--gene-list',
        help='Give the path to a .txt file with target genes.'
    )

    args = parser.parse_args()
    out_dir = args.target_junction_directory
    now = datetime.now().strftime('%m-%d-%Y_%H.%M.%S')
    batch_num = args.batch_number
    genelist = args.gene_list

    try:
        jx_df = pd.read_table(os.path.join(out_dir, _ANNOTATED_S_FILE))
    except FileNotFoundError:
        jx_df = pd.read_table(os.path.join(out_dir, _SHORTLIST_FILE))
    filter_jxs(jx_df, out_dir, batch_num, genelist)
