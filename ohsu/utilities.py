import glob
from intervaltree import IntervalTree
import os
import pandas as pd
import re
from transcript_features import Transcript, Exons, jx_modified_tx_name

_PROJECT_PATH = os.path.dirname(os.path.realpath(__file__))
_FILES_DIR = 'files'
_SCRIPTS_DIR = 'scripts'
_GDC_TCGA_INFO = os.path.join(
    _PROJECT_PATH, 'files', 'gdc_filename_mapper',
    'gdc_sample_sheet.2020-02-27.all.tsv'
)
_TCGA_METADATA = os.path.join(
    _PROJECT_PATH, 'files', 'TCGA_full_metadata_RNA.GDC_ID.tsv'
)
_TRANSLATION_PICKLE = os.path.join(
    _PROJECT_PATH, 'files', 'trans_ranges.pickle'
)
_SHORTLIST_JXS = os.path.join(
    _PROJECT_PATH, 'files', 'actual_shortlist_jxs.json'
)
_INTERMEDIATE_FILTER_DATA = 'intermediate_filter_data.json'
_SHORTLIST_FILE = 'all_sample_shortlist.tsv'
_ANNOTATED_S_FILE = 'annotated_' + _SHORTLIST_FILE
_FILT_ANN_S_FILE = 'filter_annotated_' + _SHORTLIST_FILE
_GTEX_ALL = 'all_GTEx'
_GTEX_RESTRICT = 'GTEx_brain_testis'
_GTEX_CORE = 'GTEx_core'
_GTEX_BR = 'GTEx_breast'
_GTEX_OV = 'GTEx_ovary'
_TCGA_ALL = 'all_TCGA'
_BRCA_NORM = 'BRCA_normal'
_FILTERED_FA_COUNT = 'frame_agnostic_neoepitope_count'
_FILTERED_IF_COUNT = 'in-frame_neoepitope_count'
_IF_PRE_EPS = 'prefiltered_in-frame_epitopes'
_PRE_IF_EP_COUNT = 'prefiltered_inframe_epitope_count'
_FILTERED_IF_EPS = 'in-frame_neoepitopes'
_IF_PEP = 'in-frame_peptide_sequence'
_IH_IF_PEPS = 'hanging_txs_included_inframe_pepseqs'
_IF_BIEX = 'in-frame_all-transcript_biexons'
_IF_NH_BIEX = 'in-frame_nonhanging-tx_biexons'
_MOTIF_FILTER = 'motif_filter'

_MOTIF = 'motif'
_MOTIF_OPTS = [1, 0]

_MIN_SAMP_READS = 'min_sample_reads'
_MSR_OPTS = [0]

_MIN_COHORT_SAMPS = '#_of_cohort_samples'
_MCS_OPTS = {'BRCA': [0, 1, 5], 'OV': [0, 1, 5]}

_MIN_CHRT_SAMP_READS = 'reads_per_cohort_sample'
_MCSR_OPTS = [0, 2]
_MCSR_SHORT = _MCSR_OPTS[1:]

_MAX_NORMAL_READS = 'reads_per_normal_sample'
_MNR_OPTS = [0, 1, 3]
_MAX_MNR = _MNR_OPTS[-1]

_MAX_NORMAL_SAMPS = '#_normal_samples_allowed'
_MNS_OPTS = [1, 2, 10]

_NORM_FILTER = 'normal_cohort_id'
_FILTER_DICT = {'paired': 0, 'GTEx': 1}


_NF_OPTS = _FILTER_DICT.keys()
_ETH_FILTER_DICT = {'matchedNormals': 0, 'Gtex': 1}

_FILTER_PARAMS = {
    _MOTIF_FILTER: _MOTIF_OPTS, _MIN_SAMP_READS: _MSR_OPTS,
    _MAX_NORMAL_READS: _MNR_OPTS, _MAX_NORMAL_SAMPS: _MNS_OPTS,
    _NORM_FILTER: _NF_OPTS, _MIN_COHORT_SAMPS: _MCS_OPTS,
    _MIN_CHRT_SAMP_READS: _MCSR_OPTS
}
_FILTER_PARAMS_LONG = {
    _MOTIF_FILTER: _MOTIF_OPTS, _MIN_SAMP_READS: _MSR_OPTS,
    _MAX_NORMAL_READS: [0] + _MNR_OPTS, _MAX_NORMAL_SAMPS: [0] + _MNS_OPTS,
    _NORM_FILTER: _NF_OPTS, _MIN_COHORT_SAMPS: _MCS_OPTS,
    _MIN_CHRT_SAMP_READS: _MCSR_OPTS
}
_CANON_MOTIFS = {'GTAG', 'GCAG', 'ATAC'}
_PEP_LEN = 9
_TARGET_FILE_CANCER_MAP = {
    'TCGA-AO-A0JM-01A-21R-A056-07': 'BRCA',
    'TCGA-BH-A18V-01A-11R-A12D-07': 'BRCA',
    'TCGA-A2-A0D2-01A-21R-A034-07': 'BRCA',
    'TCGA-A2-A0SX-01A-12R-A084-07': 'BRCA',
    'TCGA-C8-A12P-01A-11R-A115-07': 'BRCA',
    'TCGA-25-1319-01A-01R-1565-13': 'OV',
    'TCGA-25-1313-01A-01R-1565-13': 'OV',
    'TCGA-61-2008-01A-02R-1568-13': 'OV',
    'TCGA-24-1431-01A-01R-1566-13': 'OV',
    'TCGA-24-2298-01A-01R-1569-13': 'OV'
}
_TARGET_FILE_TO_ID_MAP = {
    '608a335d-41c4-4ba6-afe3-93a38bbc0612': 'TCGA-AO-A0JM-01A',
    '6c22b3d5-a203-43e6-8563-f78f1eba1d28': 'TCGA-BH-A18V-01A',
    'b1d17c69-5fe7-409c-b37c-a269595d4c9d': 'TCGA-BH-A18V-06A',
    '71c5ab4f-ce13-432d-9a90-807ec33cf891': 'TCGA-A2-A0D2-01A',
    '2a6f86bc-f835-40d6-8761-9e3da237c8fa': 'TCGA-A2-A0SX-01A',
    '81309f84-9f48-4626-97dc-5c5370de0cdc': 'TCGA-C8-A12P-01A',
    '4b3449af-bf5f-4da2-a9cf-5bec71f7d776': 'TCGA-25-1319-01A',
    '271d24e6-ab25-4d18-9eb8-855cecfc906d': 'TCGA-25-1313-01A',
    'b0531660-67ac-4537-82c1-dcb2586cde23': 'TCGA-61-2008-01A',
    '3dd68b83-9c11-4c4b-ba98-4fc9d8bf17b3': 'TCGA-24-1431-01A',
    '2956bd30-1b96-468f-a017-d90303098d40': 'TCGA-24-2298-01A'
}
_TARGET_ID_TO_FILE_MAP = {
    value: key for key, value in _TARGET_FILE_TO_ID_MAP.items()
}
_TARGET_SAMPLE_IDS = _TARGET_ID_TO_FILE_MAP.keys()
_TARGET_FILENAMES = _TARGET_FILE_TO_ID_MAP.keys()
# storage file names for consistency across scripts
_NORMAL_SAMPLE_ID_JSON = 'gtex_tcga_normal_sampleID_map.json'
_TUMOR_SAMPLE_ID_JSON = 'tcga_tumor_sampleID_map.json'

_TOT_JXS = 'total_jxs'
_NORMS = 'normal'
_COHORT = 'cohort'
_DF_FILE_TAGS = {
    _NORMS: 'norm', _COHORT: 'chrt', _MOTIF_FILTER: 'mot',
    _MIN_SAMP_READS: 'samp'
}
_PARAM_COLS = {
    _MAX_NORMAL_SAMPS: 'max # normal\nsamples allowed',
    _MOTIF_FILTER: 'motif\nfilter',
    _MAX_NORMAL_READS: 'max # normal\nreads allowed',
    _MIN_SAMP_READS: 'min # reads\nin sample',
    _NORM_FILTER: 'normal cohort for\nbackground filter',
    # _MIN_COHORT_READS: 'min required reads\nacross cancer cohort',
    _MIN_CHRT_SAMP_READS: 'min required reads\nper cancer cohort sample',
    _MIN_COHORT_SAMPS: 'min required\ncohort samples',
}

_ASN_COUNT_COLS = {
    _MAX_NORMAL_SAMPS: 'Filter_Sample_Cohort_CohortBackground',
    _MOTIF_FILTER: 'Filter_Motif',
    _MAX_NORMAL_READS: 'Filter_Sample_Cohort_CohortBackground',
    _MIN_SAMP_READS: 'Filter_Sample',
    _NORM_FILTER: 'Filter_Sample_Cohort_CohortBackground',
    _MIN_CHRT_SAMP_READS: 'Filter_Sample_Cohort',
    _MIN_COHORT_SAMPS: 'Filter_Sample_Cohort',
    _COHORT: 'Filter_Sample_Cohort',
    _NORMS: 'Filter_Sample_Cohort_CohortBackground'
}
_UNIPROT_FINAL_COL = 'Filter_Sample_Cohort_CohortBackground_Uniprot'
_VALIDATED_FINAL_COL = 'validated_kmer_count'
_HEATMAP_COLORBAR_NAMES = {
    _MOTIF_FILTER: 'filtered on\nsplice motifs',
    _MIN_SAMP_READS: 'filtered on\nmin sample reads',
    _COHORT: 'filtered on\ncohort foreground',
    _NORMS: 'filtered on\nnormal background'
}
_FILT_FILE_TARGET_COL = 'Filter_Sample_Cohort_CohortBackground_Uniprot'

_FILT_FILENAME_MAP = {
    'matchedNormals': 'paired', 'GTEXcore': 'core_GTEx', 'GtexTcga': 'All'
}

def info_from_SplAdder_line(SplAdder_line):
    """Returns junction (str) and coverage (int) from one SplAdder jx line"""
    items = SplAdder_line.decode('utf-8').strip().split('\t')
    chrom, strand, left, right = items[0:4]
    # Add 1 to left to match SJ.out format
    jx = ';'.join([chrom, str(int(left) + 1), right, strand])
    coverage = int(items[4])
    return jx, coverage


def under_x_count_key(tissue, cutoff_val):
    return '{}_under{}'.format(tissue, cutoff_val)


def over_x_count_key(tissue, cutoff_val):
    return '{}_over{}'.format(tissue, cutoff_val)


def filter_column(sample, parameter, filter_type):
    return '{}-{}-{}'.format(sample, parameter, filter_type)


def get_trans_id(gtf_info):
    return re.sub(r".*transcript_id \"(\S+?)\"[;].*", r"\1", gtf_info)


def normal_filter(nonzero1, over_max1, max_sample_count, nonzero2=0,
                  over_max2=0, nonzero3=0, over_max3=0):
    value = int(
        ((nonzero1 + nonzero2 + nonzero3) < max_sample_count)
        and
        ((over_max1 + over_max2 + over_max3) == 0)
    )
    return value


def txs_from_gtf(gtf_file, genelist=[]):
    """Extracts splice site annotation from .gtf file

    This function is a modified version of one that is part of HISAT.
    Copyright 2014, Daehwan Kim <infphilo@gmail.com>

    HISAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    HISAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HISAT.  If not, see <http://www.gnu.org/licenses/>.
    """
    trans = {}
    trans_by_gene = {}
    trans_by_geneid = {}
    protein_coding_txs = set()
    protein_coding_genes = set()
    pc_genes_withpc_txs = set()
    tx_cds = {}
    tx_features = {}
    gene_tx_map = {}
    geneid_tx_map = {}
    annotations = {}
    searchable_genename_tree = {}
    searchable_geneid_tree = {}
    protein_coding_both_list = []
    with open(gtf_file) as gtf:
        # Parse valid exon lines from the annotation file into a dict by
        # transcript_id
        for line in gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()

            try:
                item = line.split('\t')
                chrom, left, right, strand = item[0], item[3], item[4], item[6]
                feature, values = item[2], item[8]
            except ValueError:
                continue

            left, right = int(left), int(right)
            needed_features = {'exon', 'start_codon', 'stop_codon', 'CDS'}
            if feature not in needed_features or left >= right:
                continue

            vals_dict = dict(re.findall(r"(\w+)\s\"(\S+?)\"[;]", values))

            if 'gene_id' not in vals_dict or 'transcript_id' not in vals_dict:
                continue

            if genelist:
                if 'gene_name' not in vals_dict:
                    continue
                if vals_dict['gene_name'] not in genelist:
                    continue

            curr_tx = vals_dict['transcript_id']

            if vals_dict.get('transcript_type') == 'protein_coding':
                if vals_dict.get('gene_type') == 'protein_coding':
                    if curr_tx not in protein_coding_txs:
                        pc_genes_withpc_txs.add(vals_dict['gene_id'])
                        protein_coding_txs.add(curr_tx)
                        protein_coding_both_list.append(curr_tx)

            if vals_dict.get('gene_type') == 'protein_coding':
                protein_coding_genes.add(vals_dict['gene_name'])
            try:
                tx_features[curr_tx].add(feature)
            except KeyError:
                tx_features[curr_tx] = set(feature)

            if curr_tx not in gene_tx_map:
                if 'gene_name' in vals_dict:
                    gene_tx_map[curr_tx] = vals_dict['gene_name']
                if 'gene_id' in vals_dict :
                    geneid_tx_map[curr_tx] = vals_dict['gene_id']

            if feature == 'CDS':
                shift = int(item[7])
                try:
                    tx_cds[curr_tx][1].append([left, right, shift])
                except KeyError:
                    tx_cds[curr_tx] = [strand, [[left, right, shift]]]
            else:
                try:
                    trans[curr_tx][2].append([left, right, feature])
                except KeyError:
                    trans[curr_tx] = [chrom, strand, [[left, right, feature]]]

    print('# of protein coding transcripts in protein coding genes:')
    print(len(protein_coding_both_list))
    print('# of protein coding genes:')
    print(len(protein_coding_genes))
    print('# of protein coding genes with protein coding transcripts:')
    print(len(pc_genes_withpc_txs))
    with open('gencodev32_proteincodingtranscripts.txt', 'w') as output:
        for tx in protein_coding_both_list:
            output.write('{}\n'.format(tx))

    with open('gencodev32_proteincodinggenes.txt', 'w') as output:
        for gene in pc_genes_withpc_txs:
            output.write('{}\n'.format(gene))

    for tx, features in tx_features.items():
        if 'start_codon' not in features:
            try:
                strand = tx_cds[tx][0]
            except KeyError:
                del trans[tx]
                continue

            curr_cds = tx_cds[tx][1]
            if strand == '-':
                curr_cds.sort(key=lambda x: x[1], reverse=True)
                right = curr_cds[0][1] - curr_cds[0][2]
                left = right - 2
            else:
                curr_cds.sort(key=lambda x: x[0])
                left = curr_cds[0][0] + curr_cds[0][2]
                right = left + 2

            trans[tx][2].append([left, right, 'start_codon'])

    mod_ex = 'modified_exon'
    for tran, [chrom, strand, features] in trans.items():
        features.sort(key=lambda x: x[0])
        tmp_features = features
        for i, feature in enumerate(features):
            if feature[2] == 'start_codon' and i != 0:
                if strand == '+' and feature[1] < features[i - 1][1]:
                    curr_feature = [[feature[0], features[i-1][1], mod_ex]]
                    tmp_features = curr_feature + features[i + 1:]
                elif strand == '-' and feature[0] > features[i - 1][0]:
                    curr_feature = [[features[i - 1][0], feature[1], mod_ex]]
                    tmp_features = features[:i - 1] + curr_feature
        trans[tran] = [chrom, strand, tmp_features]

    # Sort exons and merge where separating introns are <=5 bps
    sc = 'stop_codon'
    for tran, [chrom, strand, exons] in trans.items():
        exons.sort(key=lambda x: x[0])
        tmp_exons = [exons[0]]
        for i in range(1, len(exons)):
            if exons[i][0] - tmp_exons[-1][1] <= 5:
                if exons[i][2] != sc and tmp_exons[-1][2] != sc:
                    tmp_exons[-1][1] = exons[i][1]
            else:
                tmp_exons.append(exons[i])
        trans[tran] = [chrom, strand, tmp_exons]
        gene = gene_tx_map.get(tran, 'other')
        gene_id = geneid_tx_map.get(tran, 'other')
        if tran in protein_coding_txs:
            try:
                trans_by_gene[gene][tran] = [chrom, strand, tmp_exons]
            except KeyError:
                trans_by_gene[gene] = {tran: [chrom, strand, tmp_exons]}
            try:
                trans_by_geneid[gene_id][tran] = [chrom, strand, tmp_exons]
            except KeyError:
                trans_by_geneid[gene_id] = {tran: [chrom, strand, tmp_exons]}

    for gene, transcripts in trans_by_gene.items():
        left = 10000000000000
        right = 0
        chrom = ''
        strand = ''
        for tx, info in transcripts.items():
            chrom, strand, exons = info
            if chrom not in searchable_genename_tree:
                searchable_genename_tree[chrom] = {}
            if strand not in searchable_genename_tree[chrom]:
                searchable_genename_tree[chrom][strand] = IntervalTree()
            for exon in exons:
                left = min(left, exon[0])
                right = max(right, exon[1])
        searchable_genename_tree[chrom][strand][left:right+1] = gene

    for gene, transcripts in trans_by_geneid.items():
        left = 10000000000000
        right = 0
        chrom = ''
        strand = ''
        for tx, info in transcripts.items():
            chrom, strand, exons = info
            if chrom not in searchable_geneid_tree:
                searchable_geneid_tree[chrom] = {}
            if strand not in searchable_geneid_tree[chrom]:
                searchable_geneid_tree[chrom][strand] = IntervalTree()
            for exon in exons:
                left = min(left, exon[0])
                right = max(right, exon[1])
        searchable_geneid_tree[chrom][strand][left:right+1] = gene

    # Calculate and return unique junctions
    junctions = set()
    for chrom, strand, exons in trans.values():
        for i in range(1, len(exons)):
            junctions.add((chrom, exons[i - 1][1], exons[i][0], strand))
    junctions = sorted(junctions)

    for chrom, left, right, strand in junctions:
        if chrom not in annotations:
            annotations[chrom] = {}
            annotations[chrom]['+'] = {'full': [], 'fivepr': [], 'threepr': []}
            annotations[chrom]['-'] = {'full': [], 'fivepr': [], 'threepr': []}
        annotations[chrom][strand]['full'].append(
            ';'.join([chrom, str(left), str(right), strand]))
        if strand == '+':
            five_site = str(left)
            three_site = str(right)
        else:
            five_site = str(right)
            three_site = str(left)
        if five_site not in annotations[chrom][strand]['fivepr']:
            annotations[chrom][strand]['fivepr'].append(five_site)
        if three_site not in annotations[chrom][strand]['threepr']:
            annotations[chrom][strand]['threepr'].append(three_site)

    to_return = (
        annotations, trans_by_gene, searchable_genename_tree,
        searchable_geneid_tree
    )
    return to_return


def txs_to_tree(gene_txs):
    searchable_tree = {}
    for gene_name, transcripts in gene_txs.items():
        if gene_name not in searchable_tree:
            searchable_tree[gene_name] = IntervalTree()
        for tx_id, info in transcripts.items():
            for exon in info[2]:
                left = exon[0]
                right = exon[1] + 1
                searchable_tree[gene_name][left:right] = tx_id
    return searchable_tree


def find_instream_exons(jx, genes, exon_interval_tree):
    fivepr_hit = 0
    threepr_hit = 0
    three_tx_set = set()
    five_tx_set = set()
    chrom, left, right, strand = jx.split(';')
    if strand == '+':
        fivepr = int(left) - 1
        threepr = int(right) + 1
    else:
        fivepr = int(right) + 1
        threepr = int(left)
    for gene in genes.split(','):
        try:
            overlap_five = exon_interval_tree[gene].overlaps(fivepr)
            overlap_three = exon_interval_tree[gene].overlaps(threepr)
        except KeyError:
            overlap_five = 0
            overlap_three = 0
        if overlap_five:
            fivepr_hit = 1
            for info in exon_interval_tree[gene][fivepr]:
                five_tx_set.add(info[2])
        if overlap_three:
            threepr_hit = 1
            for info in exon_interval_tree[gene][threepr]:
                three_tx_set.add(info[2])
    five_txs = ';'.join(list(five_tx_set))
    three_txs = ';'.join(list(three_tx_set))
    return fivepr_hit, five_txs, threepr_hit, three_txs


def jx_gene_overlap(junction, cds_tree):
    """Check found junctions for coding region overlap

    Input:
        junction information: chromosome, left and right boundaries of the
            junction, and strand
        CDS tree containing coding regions, created by gtf_to_cds
        dictionary mapping gene ids from .gtf file to gene names

    Checks to see whether either junction right overlaps coding regions.  If
    either or both do, collects gene ids and names for the overlapping genes.
    If both sides overlap, checks to see if any of the genes are the same.

    Returns a semicolon-joined string of gene names overlapping the 5' and 3'
    junction ends
    """
    left_names = set()
    right_names = set()
    chrom, left, right, strand = junction.split(';')
    left = int(left) - 1
    right = int(right) - 1
    try:
        jx_start = cds_tree[chrom][strand].overlaps(left)
        jx_stop = cds_tree[chrom][strand].overlaps(right)
    except KeyError:
        return ';'

    if jx_start or jx_stop:
        for start_set in list(cds_tree[chrom][strand][left]):
            left_names.add(start_set[2])
        for stop_set in list(cds_tree[chrom][strand][right]):
            right_names.add(stop_set[2])

    left_names = ','.join(list(left_names))
    right_names = ','.join(list(right_names))

    return ';'.join([left_names, right_names])


def build_novel_transcript_dictionary(gtf_file, target_txs, verbose=False):
    tx_dict = {}
    transcripts_of_interest = set(target_txs.keys())
    prev_line = ''
    prev_feature = ''
    with open(gtf_file) as open_gtf:
        for line in open_gtf:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '#' in line:
                line = line.split('#')[0].strip()
            try:
                item = line.split('\t')
                chrom, _, feature, start, stop, _, strand, read_frame, _ = item
            except ValueError:
                continue
            curr_tx = get_trans_id(line)
            if curr_tx in transcripts_of_interest:
                if feature == 'transcript' and curr_tx not in tx_dict.keys():
                    seleno = 'seleno' in line
                    tx_dict[curr_tx] = Transcript(
                        curr_tx, chrom, strand, seleno
                    )
                if feature == 'exon':
                    tx_dict[curr_tx].add_exon(start, stop, line)
                elif feature in ['start_codon', 'stop_codon', 'CDS']:
                    if feature == prev_feature:
                        prev_line = line
                        continue
                    else:
                        tx_dict[curr_tx].add_feature_to_exon(
                            feature, start, stop, read_frame
                        )
            prev_feature = feature
            if prev_line:
                item = line.split('\t')
                chrom, _, feature, start, stop, _, strand, read_frame, _ = item
                curr_tx = get_trans_id(line)
                tx_dict[curr_tx].add_feature_to_exon(
                    feature, start, stop, read_frame
                )
                prev_line = ''

    for tx, junction_list in target_txs.items():
        for jx in junction_list:
            newkey = jx_modified_tx_name(tx, jx)
            if verbose:
                print('original transcript for {}:'.format(tx))
                tx_dict[tx].print()
            if newkey in tx_dict:
                continue
            try:
                tx_dict[newkey] = tx_dict[tx].copy_transcript(newkey)
            except KeyError:
                print('key error in tx dict creation!')
                print(tx, newkey)
                continue
            if verbose:
                print('copied transcript for {}:'.format(newkey))
                tx_dict[newkey].print()
            chrom, left, right, strand = jx.split(';')
            # update junction coordinates to exon coordinates:
            left = str(int(left) - 1)
            right = str(int(right) + 1)
            jx = ';'.join([chrom, left, right, strand])
            tx_dict[newkey].insert_junction(jx, verbose=verbose)
            if verbose:
                print('junction-inserted transcript for {}:'.format(newkey))
                tx_dict[newkey].print()
    return tx_dict


def load_filter_dataframe(result_dir, sample):
    file = os.path.join(result_dir, '*filtered_df_{}*'.format(sample))
    fileglob = glob.glob(file)[0]
    if os.path.basename(fileglob).startswith('G'):
        eth = True
        separator=','
        renames = {
            'Filter_Sample_Cohort_CohortNormal_Uniprot': _FILT_FILE_TARGET_COL
        }
    else:
        eth = False
        separator='\t'
        renames = {}
    statdf = pd.read_table(fileglob, sep=separator)
    if eth:
        statdf = statdf.loc[~statdf['normal_cohort_id'].isin({
            'GtexTcga1', 'GtexTcga2', 'GtexTcga3', 'GtexTcga4', 'GtexTcga5',
            'GtexTcga6', 'GtexTcga7', 'GtexTcga8', 'GtexTcga9', 'GtexTcga0',
            'GTEX'
        })].copy()
        statdf.rename(renames, axis=1, inplace=True)
        statdf[_NORM_FILTER] = statdf[_NORM_FILTER].apply(
            lambda x: _FILT_FILENAME_MAP[x]
        )
        nonemap = {'None': 0}
        statdf[_MIN_COHORT_SAMPS] = statdf[_MIN_COHORT_SAMPS].apply(
            lambda x: nonemap.get(x, x)
        ).astype(int)
        statdf[_MIN_CHRT_SAMP_READS] = statdf[_MIN_CHRT_SAMP_READS].apply(
            lambda x: nonemap.get(x, x)
        ).astype(float)
        statdf[_MAX_NORMAL_SAMPS] = statdf.apply(
            lambda x: 0 if x[_MAX_NORMAL_READS] == 0 else x[_MAX_NORMAL_SAMPS],
            axis=1
        )
        statdf = statdf.loc[statdf[_MAX_NORMAL_SAMPS] != 1].copy()
        statdf[_MAX_NORMAL_READS] = statdf[_MAX_NORMAL_READS].astype(int)
        statdf[_MIN_SAMP_READS] = statdf[_MIN_SAMP_READS].astype(int)
        statdf[_MIN_CHRT_SAMP_READS] = statdf[_MIN_CHRT_SAMP_READS].astype(int)
    return statdf


def load_validated_dataframe(result_dir, sample):
    file = os.path.join(result_dir, '*validated_df_{}*'.format(sample))
    fileglob = glob.glob(file)[0]
    if os.path.basename(fileglob).startswith('G'):
        eth = True
        separator=','
        renames = {
            'Filter_Sample_Cohort_CohortNormal_Uniprot': _FILT_FILE_TARGET_COL
        }
    else:
        eth = False
        separator='\t'
        renames = {}
    statdf = pd.read_table(fileglob, sep=separator)
    if eth:
        statdf = statdf.loc[~statdf['normal_cohort_id'].isin({
            'GtexTcga1', 'GtexTcga2', 'GtexTcga3', 'GtexTcga4', 'GtexTcga5',
            'GtexTcga6', 'GtexTcga7', 'GtexTcga8', 'GtexTcga9', 'GtexTcga0',
            'GTEX'
        })].copy()
        statdf.rename(renames, axis=1, inplace=True)
        statdf[_NORM_FILTER] = statdf[_NORM_FILTER].apply(
            lambda x: _FILT_FILENAME_MAP[x]
        )
        nonemap = {'None': 0}
        statdf[_MIN_COHORT_SAMPS] = statdf[_MIN_COHORT_SAMPS].apply(
            lambda x: nonemap.get(x, x)
        ).astype(int)
        statdf[_MIN_CHRT_SAMP_READS] = statdf[_MIN_CHRT_SAMP_READS].apply(
            lambda x: nonemap.get(x, x)
        ).astype(float)
        statdf[_MAX_NORMAL_SAMPS] = statdf.apply(
            lambda x: 0 if x[_MAX_NORMAL_READS] == 0 else x[_MAX_NORMAL_SAMPS],
            axis=1
        )
        statdf = statdf.loc[statdf[_MAX_NORMAL_SAMPS] != 1].copy()
        statdf[_MAX_NORMAL_READS] = statdf[_MAX_NORMAL_READS].astype(int)
        statdf[_MIN_SAMP_READS] = statdf[_MIN_SAMP_READS].astype(int)
        statdf[_MIN_CHRT_SAMP_READS] = statdf[_MIN_CHRT_SAMP_READS].astype(int)
    return statdf
