#!/usr/bin/env python3

"""
collect_shortlist_jxs.py
Python 3.6 code for collecting shortlist TCGA junctions from SJ.out.tab
junction files called by the STAR aligner.

Junctions provided by Ratsch Lab, Biomedical Informatics at ETH Zurich.
Alignment parameters:
STAR --genomeDir GENOME_INDEX --readFilesIn FQ1 FQ2 --runThreadN 8
--sjdbOverhang 100 --outFileNamePrefix OUT_PREFIX
--outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20
--outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000
--sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory
--readFilesCommand zcat --outFilterMatchNminOverLread 0.33
--outFilterScoreMinOverLread 0.33 --outSAMstrandField intronMotif
--outSAMmode Full --outSAMattributes NH HI NM MD AS XS --outSAMunmapped Within
--limitSjdbInsertNsj 2000000 --outSAMtype BAM Unsorted
--outSAMheaderHD @HD VN:1.4 --twopassMode Basic --outSAMmultNmax 1

Alignment output .bam files processed with SplAdder yield the final junction
inputs; total coverage count files provided by Ratsch Lab, BMI, ETH Zurich.

GTEx phenotype files with recount2-generated phenotype information can be
obtained from https://jhubiostatistics.shinyapps.io/recount/: click on GTEx
and then the "link" in the phenotype column.

TCGA File Mapper, _GDC_TCGA_INFO, comes from repository data:
https://portal.gdc.cancer.gov/repository
provided by Ratsch Lab, BMI, ETH Zurich

Sample run: home_dir %
time python ../junction-filtering/scripts/collect_shortlist_jxs.py
-T files/tcga_spladder/ -G files/gtex_spladder/
-t files/alignment.coverage_stats_tcga.tsv.gz
-g files/alignment.coverage_stats_gtex.tsv.gz
-p ../immunotherapy/files/SRP012682.tsv

"""
import argparse
import glob
import gzip
import os
import pandas as pd
import sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from utilities import _GDC_TCGA_INFO, _TARGET_FILE_CANCER_MAP, _NF_OPTS
from utilities import _SHORTLIST_FILE, _GTEX_CORE, _GTEX_RESTRICT, _TCGA_ALL
from utilities import _GTEX_BR, _MAX_MNR, over_x_count_key, under_x_count_key
from utilities import info_from_SplAdder_line, _GTEX_OV, _MNR_OPTS
from utilities import _TCGA_METADATA, _MCSR_OPTS, _MCSR_SHORT

# Previous coverage multiplier based on AUC normalization:
# _COV_MULT = 1e10
# Current multiplier based on upper quartile normalization:
_COV_MULT = 4e5
_RESTRICTED_NORMALS = {'Brain', 'Testis'}
_CANCER_MATCH = {'Breast': 'BRCA', 'Ovary': 'OV'}
_KEY_MATCH = {'Breast': _GTEX_BR, 'Ovary': _GTEX_OV}
_NORMALS = [_GTEX_BR, _GTEX_OV, _GTEX_CORE, _GTEX_RESTRICT, _TCGA_ALL]
_NORMAL_DICT_FLAGS = {
    over_x_count_key(tis, val) for tis in _NORMALS for val in [0] + _MNR_OPTS
}
_COHORT_DICT_FLAGS = {
    '{}_over_{}'.format(cancer, val)
    for cancer in ['BRCA', 'OV'] for val in _MCSR_OPTS
}


def parse_covs(gtex_cov, tcga_cov):
    """Reads GTEx and TCGA coverage files maps IDs to total sample coverage

    Input:
        gtex_auc_covs (str): path to gtex sample coverage file from ETH
            alignment.coverage_stats_gtex.tsv.gz
        tcga_cov (str): path to tcga sample coverage file from ETH
            alignment.coverage_stats_tcga.tsv.gz
        file_to_id (dict): maps TCGA filenames to IDs

    Returns dict mapping GTEx/TCGA ID to total sample coverage
    """
    cov_map = {}
    for cov_file in [gtex_cov, tcga_cov]:
        with open(cov_file) as cov_lines:
            next(cov_lines)
            cov_index = 1
            for line in cov_lines:
                line = line.strip().split('\t')
                strid, coverage = line[0], line[cov_index]
                try:
                    cov = int(float(coverage))
                except ValueError:
                    continue
                if cov_file == gtex_cov and strid.endswith('all'):
                    strid = strid.strip('all')
                if cov_file == tcga_cov:
                    strid = '-'.join([
                        strid[0:4], strid[4:6], strid[6:10], strid[10:13],
                        strid[13:16], strid[16:20], strid[20:22]
                    ])
                cov_map[strid] = cov
    return cov_map


def parse_phenotypes(tissue_phen, cancer_phen):
    """Parses phenotype files to generate sample type and id name maps.

    Input:
        tissue_phen (str): path to recount-generated file w/ sample SRA numbers
            and tissue types for normal tissue samples
        cancer_phen (str): gdc file, stored in repo: see _GDC_TCGA_INFO above

    Parses GTEx and TCGA info files.
    Returns two dictionaries: 1) tissue_map maps file ID for normal samples to
    tissue type; 2) file_to_id maps filenames to TCGA sample IDs.
    """
    gtex_phen = pd.read_table(tissue_phen, usecols=['sampid', 'smts', 'run'])
    new_names = {
        'smts': 'sample_type', 'run': 'universal_id', 'sampid': 'gtex_id'
    }
    gtex_phen.rename(new_names, axis='columns', inplace=True)
    gtex_phen = gtex_phen[gtex_phen.gtex_id.str.startswith('GTEX')]
    gtex_phen.dropna(subset=['sample_type'], inplace=True)
    tissue_id_index = gtex_phen.set_index('universal_id')['sample_type']
    tcga_cols = [
        'study', 'tcga_id', 'is_normal'
    ]
    tcga_phen = pd.read_table(cancer_phen, usecols=tcga_cols)
    new_names = {
        'study': 'abbr', 'is_normal': 'normal'
    }
    tcga_phen.rename(new_names, axis='columns', inplace=True)
    target_samps = tcga_phen.loc[
        tcga_phen['abbr'].isin(['BRCA', 'OV'])
    ]
    # Note: column 'is_normal' is read in as a bool
    target_samps = target_samps[target_samps['normal'] == False]
    cancer_cohort_dict = target_samps.set_index('tcga_id')['abbr'].to_dict()
    tcga_phen.sort_values(by=['tcga_id'], inplace=True)
    tcga_phen = tcga_phen[tcga_phen['normal'] == True]
    cancer_id_index = tcga_phen.set_index('tcga_id')['abbr']
    tissue_map = pd.concat([tissue_id_index, cancer_id_index]).to_dict()
    return tissue_map, cancer_cohort_dict


def collect_target_jxs(tcga_dir, cov_map):
    """Collects all junctions from target BRCA and OV samples for filtering.

    Input:
    gtex_dir (str): path to directory w/ TCGA sample junction files
    auc_cov_map (dict): dict mapping TCGA IDs to total sample coverage values

    Reads junction files for target samples and creates two dicts:
        - jxs_sets contains the sets of all junctions in OV and BRCA
        - all_jx_info has all junctions as keys,  with coverage for each sample

    Returns tuple of (jxs_sets, all_jx_info)
    """
    jx_sets = {'OV': set(), 'BRCA': set()}
    all_jx_info = {}
    for f_id, cancer in _TARGET_FILE_CANCER_MAP.items():
        try:
            file = glob.glob(os.path.join(tcga_dir, '{}*'.format(f_id)))[0]
        except IndexError:
            print('{} file not present, continuing'.format(f_id))
            continue
        with gzip.open(file) as sample_jxs:
            next(sample_jxs)
            for line in sample_jxs:
                # COUNTING UNIQUE-MAPPER COVERAGE ONLY
                jx, cov = info_from_SplAdder_line(line)
                jx_sets[cancer].add(jx)
                try:
                    all_jx_info[jx][f_id] = _COV_MULT * cov / cov_map[f_id]
                except KeyError:
                    all_jx_info[jx] = {}
                    all_jx_info[jx][f_id] = _COV_MULT * cov / cov_map[f_id]
    return (jx_sets, all_jx_info)


def cancer_cohort_jxs(target_jx_info, cohort_sample_map, tcga_dir, cov_map):
    jx_cancer_dict = {}
    target_jxs = target_jx_info[0]['OV'].union(target_jx_info[0]['BRCA'])
    for f_id, cancer in cohort_sample_map.items():
        try:
            file = glob.glob(os.path.join(tcga_dir, '{}*'.format(f_id)))[0]
        except IndexError:
            print('{} file not present, continuing'.format(f_id))
            continue
        with gzip.open(file) as sample_jxs:
            next(sample_jxs)
            for line in sample_jxs:
                # COUNTING UNIQUE-MAPPER COVERAGE ONLY
                jx, cov = info_from_SplAdder_line(line)
                if jx not in target_jxs:
                    continue
                if jx not in jx_cancer_dict.keys():
                    # Begin at -1 to remove target sample from count
                    jx_cancer_dict[jx] = {
                        key: 0 for key in _COHORT_DICT_FLAGS
                    }
                cov = _COV_MULT * cov / cov_map[f_id]
                jx_cancer_dict[jx]['{}_over_0'.format(cancer)] += 1
                for cov_cutoff in _MCSR_SHORT:
                    if cov >= cov_cutoff:
                        key = '{}_over_{}'.format(cancer, cov_cutoff)
                        jx_cancer_dict[jx][key] += 1
    return jx_cancer_dict


def junctions_from_normal_samples(files, tissue_type_map, cov_map, jx_sets,
                                  tissue_target_samps):
    """Collects junction info from normal samples for cancer junction filtering

    Input:
    files (list): list of paths to all GTEx and TCGA files
    target_ids (dict): maps IDs to tissue types, from parse_phenotypes
    cov_map (dict): maps IDs to total sample coverage, from parse_gtex_auc_covs
    jx_sets(dict): dict containing cancer:set of all junctions for OV & BRCA

    Reads junctions from each normal sample file in GTEx and TCGA directories.
    Junction data is summarized for future junction filtering.

    Returns dict of junctions w/ aggregated coverage cutoff sample count data
    """
    all_target_jxs = jx_sets['OV'].union(jx_sets['BRCA'])
    normal_jx_info = {}
    for filename in files:
        print('\nnew file:')
        id = os.path.basename(filename).split('.')[0]
        try:
            tissue = tissue_type_map[id]
            samp_cov = cov_map[id]
            print('{}: {}'.format(id, tissue))
        except KeyError:
            print('not a normal sample ({}) continuing'.format(id))
            continue
        if id.startswith('SRR'):
            if id not in tissue_target_samps:
                continue
            if tissue in _RESTRICTED_NORMALS:
                proj_flag = _GTEX_RESTRICT
            else:
                proj_flag = _GTEX_CORE
        else:
            proj_flag = _TCGA_ALL

        with gzip.open(filename) as sample_jxs:
            next(sample_jxs)
            for line in sample_jxs:
                jx, jx_cov = info_from_SplAdder_line(line)
                if jx not in all_target_jxs:
                    continue
                if jx not in normal_jx_info:
                    normal_jx_info[jx] = {
                        header: 0 for header in _NORMAL_DICT_FLAGS
                    }
                try:
                    # We only care about breast junction coverage if the jx is
                    # also in BRCA, and ovary coverage if the jx is in OV
                    cancer = _CANCER_MATCH[tissue]
                    match = jx in jx_sets[cancer]
                    tis_key = _KEY_MATCH[tissue]
                except KeyError:
                    match = False
                coverage = _COV_MULT * jx_cov / samp_cov
                # for cov_cutoff in [0] + _MNR_OPTS:
                for cov_cutoff in _MNR_OPTS:
                    if cov_cutoff == 0:
                        jx_meets_cov_cutoff = coverage > cov_cutoff
                    else:
                        jx_meets_cov_cutoff = coverage >= cov_cutoff
                    if jx_meets_cov_cutoff:
                        dict_key = over_x_count_key(proj_flag, cov_cutoff)
                        # count all samples at each coverage cutoff level
                        normal_jx_info[jx][dict_key] += 1
                        if match:
                            # count tissue-matched samples at coverage cutoff
                            dict_key = over_x_count_key(tis_key, cov_cutoff)
                            normal_jx_info[jx][dict_key] += 1
    return normal_jx_info


def print_normal_sample_list(tissue_type_map, tcga_dir, gtex_dir):
    gtex_files = glob.glob(os.path.join(gtex_dir, '*.all*tsv.gz'))
    tcga_files = glob.glob(os.path.join(tcga_dir, '*.all*tsv.gz'))
    gtex_normal = []
    tcga_normal = []
    for filename in gtex_files:
        id = os.path.basename(filename).split('.')[0]
        try:
            _ = tissue_type_map[id]
        except KeyError:
            continue
        else:
            gtex_normal.append(id)
    with open('GTEx_sample_IDs_08-2021.csv', 'w') as output:
        for sample in gtex_normal:
            output.write('{}\n'.format(sample))
    for filename in tcga_files:
        id = os.path.basename(filename).split('.')[0]
        try:
            _ = tissue_type_map[id]
        except KeyError:
            continue
        else:
            tcga_normal.append(id)
    with open('TCGA_normal_sample_IDs_08-2021.csv', 'w') as output:
        for sample in tcga_normal:
            output.write('{}\n'.format(sample))
    return


def collect_normal_jxs(jx_info, tissue_type_map, cov_map, tcga_dir, gtex_dir,
                       cohort_jx_dict, out_path, gtex_target_samps):
    """Collects junction info from normal samples for cancer junction filtering

    Input:
    jx_info (tup): tuple of two junction info dicts, from collect_scaled_covs
    target_ids (dict): maps IDs to tissue types, from parse_phenotypes
    cov_map (dict): maps IDs to total sample coverage, from parse_gtex_auc_covs
    gtex_dir (str): path to directory w/ TCGA sample junction files
    gtex_dir (str): path to directory w/ GTEx sample junction files
    out_path (str): path to directory to store output files

    Reads junctions from each normal sample file in GTEx and TCGA directories.
    Junction data is summarized for future junction filtering.

    Returns None
    """
    cohort_jx_df = pd.DataFrame.from_dict(
        cohort_jx_dict, orient='Index'
    ).reset_index()
    cohort_jx_df.rename({'index': 'jx'}, axis=1, inplace=True)
    jx_sets, all_jx_info = jx_info
    # jx_sets: all junctions for OV, all for BRCA
    # all_jx_info: motif and ID coverage info.
    gtex_files = glob.glob(os.path.join(gtex_dir, '*.all*tsv.gz'))
    normal_jx_dict = junctions_from_normal_samples(
        # gtex_files + tcga_files, tissue_type_map, cov_map, jx_sets
        gtex_files, tissue_type_map, cov_map, jx_sets, gtex_target_samps

    )
    normal_jx_df = pd.DataFrame(normal_jx_dict).fillna(0)
    normal_jx_df = normal_jx_df.transpose().reset_index()
    normal_jx_df.rename({'index': 'jx'}, axis='columns', inplace=True)
    tumor_df = pd.DataFrame(all_jx_info).fillna(0)
    tumor_df = tumor_df.transpose().reset_index()
    tumor_df.rename({'index': 'jx'}, axis='columns', inplace=True)
    full_df = tumor_df.merge(normal_jx_df, how='left', on='jx')
    full_df = full_df.merge(cohort_jx_df, how='left', on='jx')
    full_df = full_df.fillna(0)
    with open(os.path.join(out_path, _SHORTLIST_FILE), 'w') as output:
        full_df.to_csv(output, sep='\t', index=False)
    return


def parse_gtex_samps(gtex_samp_file):
    if gtex_samp_file is None:
        return set()
    gtex_samps = pd.read_csv(gtex_samp_file, header=None, names=['samp'])
    gtex_samps['final_names'] = gtex_samps['samp'].apply(
        lambda x: x.strip('all')
    )
    return set(gtex_samps['final_names'].tolist())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collects target junctions from GTEx & TCGA SplAdder '
                    'output files.'
    )
    parser.add_argument(
        '--output-path', '-o', default='./',
        help='Path to store sample junction output.'
    )
    parser.add_argument(
        '--TCGA-junction-directory', '-T',
        help='Directory containing SplAdder files with junctions extracted '
             'from a STAR alignment run of TCGA data.'
    )
    parser.add_argument(
        '--GTEx-junction-directory', '-G',
        help='Directory containing SplAdder files with junctions extracted '
             'from a STAR alignment run of GTEx data.'
    )
    parser.add_argument(
        '--TCGA-coverage-counts', '-t',
        help='File containing total sample coverage per TCGA sample.'
    )
    parser.add_argument(
        '--GTEx-coverage-counts', '-g',
        help='File containing total sample coverage per GTEx sample.'
    )
    parser.add_argument(
        '--GTEx-phenotype-file', '-p',
        help='File with recount2-generated GTEx phenotype information via '
             'https://jhubiostatistics.shinyapps.io/recount/ > GTEx > '
             'phenotype column; click on "link".'
    )
    parser.add_argument(
        '--GTEx-target-samples', '-s',
        help='File list of target GTEx samples shared between ETH and OHSU'
    )

    args = parser.parse_args()
    out_path = args.output_path
    tcga_dir = args.TCGA_junction_directory
    gtex_dir = args.GTEx_junction_directory
    tcga_cov = args.TCGA_coverage_counts
    gtex_cov = args.GTEx_coverage_counts
    gtex_phen = args.GTEx_phenotype_file
    gtex_samp_file = args.GTEx_target_samples

    os.makedirs(out_path, exist_ok=True)

    target_gtex_samps = parse_gtex_samps(gtex_samp_file)
    print(f'{len(target_gtex_samps)} joint gtex samples')

    # Collects tissue type for GTEx and non-tumor (only) TCGA samples
    tissue_map, brca_ov_dict = parse_phenotypes(gtex_phen, _TCGA_METADATA)
    print_normal_sample_list(tissue_map, tcga_dir, gtex_dir)
    cov_map = parse_covs(gtex_cov, tcga_cov)
    tumorsample_jx_info = collect_target_jxs(tcga_dir, cov_map)
    cohort_jx_dict = cancer_cohort_jxs(
        tumorsample_jx_info, brca_ov_dict, tcga_dir, cov_map
    )
    collect_normal_jxs(
        tumorsample_jx_info, tissue_map, cov_map, tcga_dir, gtex_dir,
        cohort_jx_dict, out_path, target_gtex_samps
    )
