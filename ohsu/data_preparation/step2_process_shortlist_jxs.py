#!/usr/bin/env python3

"""
process_shortlist_jxs.py
Python 3.6 code for processing collected shortlist TCGA junctions

SAMPLE RUN: time python ../junction-filtering/scripts/process_shortlist_jxs.py
-j new_annotation_mode -g ../retained_introns/files/gencode.v34.annotation.gtf
-f ../immunotherapy/files/GRCh38.primary_assembly.genome.fa
-u files/uniprot-proteome_UP000005640.fasta


"""
import argparse
import glob
import gzip
from math import ceil, floor
import os
import pandas as pd
import pickle
import re
import subprocess as sp
import sys
sys.path.append(
    os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
)
from utilities import _SHORTLIST_FILE, _ANNOTATED_S_FILE, _TRANSLATION_PICKLE
from utilities import _FILTERED_IF_COUNT, txs_from_gtf
from utilities import txs_to_tree, find_instream_exons, jx_gene_overlap
from utilities import build_novel_transcript_dictionary, _PEP_LEN, _IF_PEP
from utilities import _IH_IF_PEPS, _IF_BIEX, info_from_SplAdder_line
from utilities import _IF_NH_BIEX, _TARGET_FILE_TO_ID_MAP
from utilities import _PRE_IF_EP_COUNT, _FILTERED_IF_EPS, _IF_PRE_EPS
from transcript_features import Transcript, jx_modified_tx_name

_ACCEPTABLE_CHROMS = {
    'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9',
    'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17',
    'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrM', 'chrX', 'chrY'
}
_REV_COMP = str.maketrans("ATCG", "TAGC")

_CODING_REGIONS = "5';3'_coding_regions"
_GENE_IDS = 'gene_ids'
_MOTIF = 'motif'
_UP_TXS = 'upstream_txs'
_DN_TXS = 'downstream_txs'
_FA_PEPS = 'frame_agnostic_pep_seqs'
_FA_BIEX = 'frame-agnostic_all-transcript_biexons'
_FA_FCOUNT = 'frame_agnostic_frame_count'
_FA_PRE_EPS = 'prefiltered_frame_agnostic_epitopes'
_PRE_FA_EP_COUNT = 'prefiltered_fa_epitope_count'
_FILTERED_FA_EPS = 'frame_agnostic_neoepitopes'
_MOD_UP_TX = 'modified_upstream_txs'
_MOD_DN_TX = 'modified_downstream_txs'
_IH_FA_PEPS = 'hanging_txs_included_all-frame_pepseqs'


def get_trans_type(gtf_info):
    return re.sub(r".*transcript_type \"(\S+?)\"[;].*", r"\1", gtf_info)


def get_gene_name(gtf_info):
    return re.sub(r".*gene_name \"(\S+?)\"[;].*", r"\1", gtf_info)


def get_gene_id(gtf_info):
    return re.sub(r".*gene_id \"(\S+?)\"[;].*", r"\1", gtf_info)


def check_annotations(junction, junction_dict):
    """Adds annotated splice junctions from .gtf file to the junction list.

    Junction column key:
        0 = neither junction side is annotated
        1 = one junction side is annotated
        2 = both junction sides are annotated, but not together
        3 = junction is fully annotated

    NOTE: right right of junction is 2 locations lower than what is annotated in
          GENCODE. To accommodate this, all junction right ends have 2 added
          before being checked against the junction annotation dictionary.

    """
    chrom, left, right, strand = junction.split(';')
    left = str(int(left) - 1)
    right = str(int(right) + 1)
    junction = chrom + ';' + left + ';' + right + ';' + strand
    try:
        if junction in junction_dict[chrom][strand]['full']:
            annotated_col = 3
        else:
            annotated_col = 0
            if strand == '+':
                five_site = left
                three_site = right
            else:
                five_site = right
                three_site = left
            if five_site in junction_dict[chrom][strand]['fivepr']:
                annotated_col += 1
            if three_site in junction_dict[chrom][strand]['threepr']:
                annotated_col += 1
    except KeyError:
        annotated_col = 0
    return annotated_col


def jx_to_motif(jx, reference_fasta, samtools_path, hyphen=False):
    """Given an input junction and reference genome, transcribes RNA sequence.

    Input:
    jx: a junction in 'chr_;left;right;strand' format, with 0-based fully
        closed coordinates, (string)
    reference_fasta: reference genome fasta file, previously sorted and
        indexed by samtools (path to fasta file, string)
    samtools_path: to be called from subprocess to collect the sequence (path
        to samtools executable, string)

    Returns a left-to-right nucleotide sequence on either side of the aberrant
    junction sufficiently long to generate the desired protein sequence.
    """
    chrom, left, right, strand = jx.split(';')
    if chrom not in _ACCEPTABLE_CHROMS:
        return ''
    left_start = left
    left_stop = str(int(left) + 1)
    left_range = chrom + ':' + left_start + '-' + left_stop
    left_output = sp.check_output(
        ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
         '{}'.format(left_range)]
    )
    left_seq = ''.join(left_output.decode("utf-8").splitlines()[1:])
    right_start = str(int(right) - 1)
    right_stop = right
    right_range = chrom + ':' + right_start + '-' + right_stop
    right_output = sp.check_output(
        ['{}'.format(samtools_path), 'faidx', '{}'.format(reference_fasta),
         '{}'.format(right_range)]
    )
    right_seq = ''.join(right_output.decode("utf-8").splitlines()[1:])
    sequence = left_seq + right_seq
    if strand == '-':
        sequence = sequence.translate(_REV_COMP)[::-1]
    if hyphen:
        sequence = sequence[:2] + '-' + sequence[-2:]
    return sequence


def kmerize_peptides(peptide_string):
    kmers = []
    try:
        peptides = peptide_string.split(';')
    except AttributeError:
        return ''

    for peptide in peptides:
        pep_len = len(peptide)
        all_kmers = [peptide[x:x + 9] for x in range(pep_len - _PEP_LEN + 1)]
        kmers.extend([
            kmer for kmer in all_kmers if 'X' not in kmer
        ])
    return kmers


def filter_epitopes(epitope_string, uniprot_set):
    ep_pass_list = []
    try:
        epitopes = epitope_string.split(';')
    except AttributeError:
        return ''
    for epitope in epitopes:
        if epitope.replace('I', 'L') not in uniprot_set:
            ep_pass_list.append(epitope)
    return ';'.join(ep_pass_list)


def parse_uniprot_fasta(uniprot_file):
    reference_proteome = set()
    kmerized_proteome = set()
    with open(uniprot_file) as uniprot:
        curr_protein = ''
        for line in uniprot:
            if line.startswith('>'):
                reference_proteome.add(curr_protein)
                curr_protein = ''
                continue
            curr_protein += line.strip().replace('I', 'L')
        reference_proteome.add(curr_protein)

    for protein in reference_proteome:
        pep_len = len(protein)
        kmerized_proteome.update(
            [protein[x:x + 9] for x in range(pep_len - _PEP_LEN + 1)]
        )
    return kmerized_proteome


def collect_target_jx_set(tcga_dir):
    """Collects all junctions from target BRCA and OV samples for filtering.

    Input:
    gtex_dir (str): path to directory w/ TCGA sample junction files
    auc_cov_map (dict): dict mapping TCGA IDs to total sample coverage values

    Reads junction files for target samples and creates two dicts:
        - jxs_sets contains the sets of all junctions in OV and BRCA
        - all_jx_info has all junctions as keys,  with coverage for each sample

    Returns tuple of (jxs_sets, all_jx_info)
    """
    jx_set = set()
    for f_id in _TARGET_FILE_TO_ID_MAP.values():
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
                jx_set.add(jx)
    return list(jx_set)


def annotate_jxs(jx_df, out_dir, ref_fasta, samtools, gtf_file, uniprot_file,
                 batch=None):
    """

    :param jx_df:o
    :param out_dir:
    :param ref_fasta:
    :param samtools:
    :param gtf_file:
    :param uniprot_file:
    :return:
    """
    # Set up data for future use
    annotations, gene_txs, genename_tree, geneid_tree = txs_from_gtf(gtf_file)
    exon_interval_tree = txs_to_tree(gene_txs)
    uniprot_proteome = parse_uniprot_fasta(uniprot_file)

    # Collect basic annotation info
    if _MOTIF not in jx_df.columns.values.tolist():
        jx_df[_MOTIF] = jx_df.jx.apply(
            lambda x: jx_to_motif(x, ref_fasta, samtools)
        )
    print('checking annotation')
    jx_df['annotation'] = jx_df.jx.apply(
        lambda x: check_annotations(x, annotations)
    )
    print('checking coding regions')
    jx_df[_CODING_REGIONS] = jx_df.jx.apply(
        lambda x: jx_gene_overlap(x, genename_tree)
    )
    print('checking CDS')
    # Determine junctions likely to be translated; collect genes & transcripts
    jx_df['has_cds'] = jx_df[_CODING_REGIONS].apply(
        lambda x:
        1 if ((x.split(';')[0] == x.split(';')[1]) and (x != ';')) else 0
    )
    full_df = jx_df
    # TODO: update to "inframe translation" and "frame agnostic translation"
    #  - for inframe - all junctions must start in exon
    #  - for all - junctions start between start and stop codons
    #     (plus: chr22 only confirm).
    print(len(jx_df))
    jx_df = jx_df.loc[jx_df['has_cds'] == 1].copy()
    print(len(jx_df))
    print('checking gene')
    jx_df[_GENE_IDS] = jx_df.jx.apply(
        lambda x: jx_gene_overlap(x, geneid_tree)
    )
    jx_df['gene'] = jx_df[_CODING_REGIONS].apply(lambda x: x.split(';')[0])
    jx_df['gene_id'] = jx_df[_GENE_IDS].apply(lambda x: x.split(';')[0])
    jx_df.drop([_GENE_IDS], inplace=True, axis=1)
    print('upstr info')
    jx_df['upstr_info'] = jx_df.apply(
        lambda x:
        find_instream_exons(x['jx'], x['gene'], exon_interval_tree),
        axis=1
    )
    jx_df['upstream_exon'] = jx_df.upstr_info.apply(lambda x: x[0])
    jx_df[_UP_TXS] = jx_df.upstr_info.apply(lambda x: x[1])
    jx_df['downstream_exon'] = jx_df.upstr_info.apply(lambda x: x[2])
    jx_df[_DN_TXS] = jx_df.upstr_info.apply(lambda x: x[3])
    jx_df.drop(['upstr_info'], inplace=True, axis=1)
    semi_df = jx_df
    print(len(jx_df))
    jx_df = jx_df.loc[jx_df['upstream_exon'] == 1].copy()
    print(len(jx_df))
    # Execute all-frame peptide translation
    target_txs = {}
    print('getting target txs')
    for index, row in jx_df.iterrows():
        for tx in row[_UP_TXS].split(';'):
            try:
                target_txs[tx].append(row['jx'])
            except KeyError:
                target_txs[tx] = [row['jx']]
    print('opening pickle')
    if os.path.exists(_TRANSLATION_PICKLE):
        with open(_TRANSLATION_PICKLE, "rb") as recover:
            Transcript.FASTA_MAP.update(pickle.load(recover))
    print('building novel tx dict')
    novel_tx_dict = build_novel_transcript_dictionary(gtf_file, target_txs)
    print('getting modified transcripts')
    jx_df[_MOD_UP_TX] = jx_df.apply(
        lambda x: [
            jx_modified_tx_name(tx, x['jx']) for tx in x[_UP_TXS].split(';')
        ], axis=1
    )
    print('DOING TRANSLATIONS OF {} JUNCTIONS'.format(len(jx_df)))
    jx_df['translation_info'] = jx_df.apply(
        lambda x: translate_jxs(
            x[_MOD_UP_TX], novel_tx_dict, ref_fasta, samtools
        ), axis=1
    )
    print('translation done, writing pickle')
    with open(_TRANSLATION_PICKLE, "wb") as output:
        pickle.dump(Transcript.FASTA_MAP, output)

    if batch is not None:
        outfile = os.path.join(
            out_dir, 'len_vs_cleavage_{}.pickle'.format(batch)
        )
    else:
        outfile = os.path.join(out_dir, 'len_vs_cleavage.pickle')
    with open(outfile, 'wb') as output:
        pickle.dump(Transcript.TRYPSIN_SITES, output)

    print('processing in-frame peptides')
    # Process in-frame peptides
    jx_df[_IF_BIEX] = jx_df.translation_info.apply(lambda x: x[3])
    jx_df[_IF_NH_BIEX] = jx_df.translation_info.apply(lambda x: x[2])
    jx_df[_IF_PEP] = jx_df.translation_info.apply(lambda x: x[0])
    jx_df[_IH_IF_PEPS] = jx_df.translation_info.apply(lambda x: x[1])
    jx_df[_IF_PRE_EPS] = jx_df[_IH_IF_PEPS].apply(
        lambda x: ';'.join(kmerize_peptides(x)) if x else ''
    )
    jx_df[_PRE_IF_EP_COUNT] = jx_df[_IF_PRE_EPS].apply(
        lambda x: len(x.split(';')) if x else 0
    )
    jx_df[_FILTERED_IF_EPS] = jx_df[_IF_PRE_EPS].apply(
        lambda x: filter_epitopes(x, uniprot_proteome) if x else ''
    )
    jx_df[_FILTERED_IF_COUNT] = jx_df[_FILTERED_IF_EPS].apply(
        lambda x: len(x.split(';')) if x else 0
    )

    print('processing frame agnostic peptides')
    # TODO: translate only chr22 frame-agnostically
    # # Process frame-agnostic peptides
    # jx_df[_FA_PEPS] = jx_df.translation_info.apply(lambda x: x[0])
    jx_df[_FA_BIEX] = jx_df.translation_info.apply(lambda x: x[4])
    jx_df['fasta_info'] = jx_df.translation_info.apply(lambda x: x[5])
    fasta_info_dict = jx_df.set_index('jx')['fasta_info'].to_dict()
    if batch is not None:
        outfile = os.path.join(
            out_dir, 'fasta_info_dict{}.pickle'.format(batch)
        )
    else:
        outfile = os.path.join(out_dir, 'fasta_info_dict.pickle')
    with open(outfile, 'wb') as output:
        pickle.dump(fasta_info_dict, output)
    jx_df.drop(['fasta_info'], axis=1, inplace=True)

    print('first merge')
    # Join all dataframes and write output
    merge_cols = semi_df.columns.values.tolist()
    jx_df = pd.merge(jx_df, semi_df, on=merge_cols, how='outer').fillna('')
    print('second merge')
    merge_cols = full_df.columns.values.tolist()
    jx_df = pd.merge(jx_df, full_df, on=merge_cols, how='outer').fillna('')
    # jx_df.drop([_FA_PRE_EPS], axis=1, inplace=True)
    jx_df.drop(['translation_info'], axis=1, inplace=True)
    print('writing file')
    if batch is not None:
        outfile = os.path.join(
            out_dir, 'batch_{}_{}'.format(batch, _ANNOTATED_S_FILE)
        )
        if batch == 0:
            print_header = True
        else:
            # To make reassembling the full file easier later
            print_header = False
    else:
        outfile = os.path.join(out_dir, _ANNOTATED_S_FILE)
        print_header = True
    with open(outfile, 'w') as output:
        jx_df.to_csv(output, sep='\t', index=False, header=print_header)
    return jx_df


def translate_jxs(target_txs, novel_tx_dict, reference_fasta, samtools_path):
    if_nonhang_biexons = set()
    all_if_biexons = set()
    if_nonhang_kmers = set()
    all_if_kmers = set()
    all_biexons = set()
    fasta_info = {}
    for tx in target_txs:
        tx_object = novel_tx_dict[tx]
        peplist, flags = tx_object.translate_jx_peptides(
            reference_fasta, samtools_path, kmer_len=_PEP_LEN
        )
        if_biex, fa_biexons, if_kmer = peplist
        if_nonhang_biexons.add(if_biex)
        fasta_info[if_biex] = {
            'jx_pos': flags[0], 'between_codons': flags[3],
            'fiveprime_end': flags[1], 'threeprime_end': flags[2]
        }
        if_nonhang_kmers.add(if_kmer)
        all_biexons.update([seq for seq in fa_biexons if len(seq) > 0])
        peplist, flags = tx_object.translate_jx_peptides(
            reference_fasta, samtools_path, kmer_len=_PEP_LEN,
            hanging_only=True
        )
        h_if_biex, h_fa_biex, if_h_kmer = peplist
        # TODO: will the fasta info potentially change between hanging/NH?
        fasta_info[h_if_biex] = {
            'jx_pos': flags[0], 'between_codons': flags[3],
            'fiveprime_end': flags[1], 'threeprime_end': flags[2]
        }
        all_if_biexons.add(h_if_biex)
        all_if_kmers.add(if_h_kmer)
        all_biexons.update([seq for seq in h_fa_biex if len(seq) > 0])
    all_if_kmers.update(if_nonhang_kmers)
    all_if_biexons.update(if_nonhang_biexons)
    to_return = (
        # In-frame, non-hanging kmers; primary peptides
        ';'.join(list(if_nonhang_kmers)).strip(';'),
        # In-frame kmers including hanging transcripts; to compare with ETH
        ';'.join(list(all_if_kmers)).strip(';'),
        # In-frame non-hanging bi-exons, for CPTAC query
        ';'.join(list(if_nonhang_biexons)).strip(';'),
        # All in-frame biexons, including hanging transcripts
        ';'.join(list(all_if_biexons)).strip(';'),
        # All-frame all-transcript biexons
        ';'.join(list(all_biexons)).strip(';'),
        fasta_info
    )
    return to_return


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Collects target junctions from GTEx & TCGA SJ.out files.'
    )
    parser.add_argument(
        '--target-junction-directory', '-j', required=True,
        help='Give the directory with the "all_sample_shortlist.tsv" file '
             'generated by collect_shortlist_jxs.py, containing collected '
             'junctions from target samples w/ aggregated normal sample '
             'counts at varying coverage cutoffs.'
    )
    parser.add_argument(
        '--gtf-file', '-g',
        help='.gtf file containing CDS annotation required to determine '
             'whether junctions are annotated and whether they occur in '
             'protein coding regions of a gene.'
    )
    parser.add_argument(
        '--reference-genome', '-f',
        help='.fa file containing the reference genome sequence. NOTE: This '
             'fasta file must previously have been indexed by running '
             '"samtools faidx <ref.fasta>".'
    )
    parser.add_argument(
        '--uniprot-proteome', '-u',
        help='.fasta file containing the uniprot reference proteome against '
             'which to query potential splicing neoepitopes. This was '
             'obtained from https://www.uniprot.org/proteomes on 9/28/2020.'
    )
    parser.add_argument(
        '--samtools-path', default='samtools',
        help='Give the path to access samtools.'
    )
    parser.add_argument(
        '--TCGA-junction-directory', '-T',
        help='Directory containing SplAdder files with junctions extracted '
             'from a STAR alignment run of TCGA data.'
    )
    parser.add_argument(
        '--batch-number', type=int,
        help='For running in batches on slurm: the slurm array task id to '
             'process correct chunk of the dataframe.'
    )
    parser.add_argument(
        '--total-batches', type=int,
        help='total number of batches for batch job on slurm'
    )

    args = parser.parse_args()
    out_dir = args.target_junction_directory
    gtf_path = args.gtf_file
    ref_fasta = args.reference_genome
    samtools = args.samtools_path
    ref_proteome = args.uniprot_proteome
    tcga_dir = args.TCGA_junction_directory
    batch_num = args.batch_number
    tot_batches = args.total_batches

    outfile = os.path.join(
        out_dir, 'batch_{}_{}'.format(batch_num, _ANNOTATED_S_FILE)
    )
    if os.path.isfile(outfile):
        print('output annotation file {} exists, exiting'.format(outfile))
        exit()

    try:
        jx_df = pd.read_table(os.path.join(out_dir, _SHORTLIST_FILE))
        if batch_num is not None:
            total_jxs = len(jx_df)
            chunksize = ceil(total_jxs / tot_batches)
            list_df = [
                jx_df[i:i + chunksize] for i in range(0, total_jxs, chunksize)
            ]
            jx_df = list_df[batch_num]
            print(
                'Batch {}: {} total jxs, {} in current batch'.format(
                    batch_num, total_jxs, len(jx_df)
                )
            )
    except FileNotFoundError:
        # Collect set of target junctions
        target_jxs = collect_target_jx_set(tcga_dir)
        target_jxs.sort()
        if batch_num is not None:
            mini_list = target_jxs[batch_num::tot_batches]
            jx_df = pd.DataFrame({'jx': mini_list})
            print(
                'Batch {}: {} total jxs, {} ({}) in current batch'.format(
                    batch_num, len(target_jxs), len(mini_list), len(jx_df)
                )
            )
        else:
            jx_df = pd.DataFrame({'jx': target_jxs})

    annotate_jxs(
        jx_df, out_dir, ref_fasta, samtools, gtf_path, ref_proteome,
        batch=batch_num
    )
