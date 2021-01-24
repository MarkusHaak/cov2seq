import os
import sys
import numpy as np
import pandas as pd
from augur.align import run as augur_align
from Bio import AlignIO, SeqIO
import edlib

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

class alignArgs:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

def mafft(fasta_fns, reference_fn, tmp_fn="clade_assignment_tmp_alignment.fasta", delete_tmp_files=True):
    if type(fasta_fns) != list:
        fasta_fns = [fasta_fns]
    aln_args = alignArgs(sequences=fasta_fns, output=tmp_fn, method='mafft',
                         reference_name=None, reference_sequence=reference_fn,
                         nthreads=4, remove_reference=False,
                         existing_alignment=False, debug=False, fill_gaps=False)
    augur_align(aln_args)
    alignment = AlignIO.read(tmp_fn, 'fasta')
    if delete_tmp_files:
        os.remove(tmp_fn)
        os.remove(tmp_fn + '.log')
        if os.path.exists(tmp_fn + ".insertions.csv"):
            os.remove(tmp_fn + ".insertions.csv")
    return alignment

def compare_consensus_to_reference(fasta_fn, reference_fn):
    target_aligned, query_aligned = align_edlib(fasta_fn, reference_fn, k=-1)
    return parse_global_alignment(target_aligned, query_aligned)

def align_edlib(fasta_fn, reference_fn, k=-1):
    '''Performs a global alignment of a single sample against the reference
       using edlib and returns the aligned target and query sequences.'''
    sample_record = next(SeqIO.parse(fasta_fn, 'fasta'))
    reference_record = next(SeqIO.parse(reference_fn, 'fasta'))
    query_seq = str(sample_record.seq).upper()
    target_seq = str(reference_record.seq).upper()
    alignment = edlib.align(query_seq, target_seq, 
                            mode='NW', task='path', k=k,
                            additionalEqualities=[('N','A'),('N','T'),('N','G'),('N','C')])
    nice = edlib.getNiceAlignment(alignment, query_seq, target_seq)
    return nice['target_aligned'], nice['query_aligned']

def parse_global_alignment(target_aligned, query_aligned):
    '''Parses the equally long strings that result of a global alignment and 
       returns a pandas DataFrame of variants, of N-masked regions as well as
       terminal alignment gaps.
       It is assumed that all characters are uppercase and that the
       alignment was produced with N matching all four canonical bases.'''
    if len(target_aligned) != len(query_aligned):
        logger.error('Sequence seems to be improperly aligned to reference since ' +\
                     'their alignment string lengths do not match')
    alignment_length = len(target_aligned)
    sites = np.empty(shape=(2, alignment_length), dtype=np.int32)
    match = np.empty(shape=(alignment_length,), dtype=np.int32)
    masked = np.empty(shape=(alignment_length,), dtype=np.int32)
    for i in range(alignment_length):
        if query_aligned[i] == target_aligned[i]:
            match[i] = 1
            masked[i] = 0
        elif query_aligned[i] == 'N' and target_aligned[i] != '-':
            match[i] = 1
            masked[i] = 1
        else:
            match[i] = 0
            masked[i] = 0

        for j,aligned in enumerate([target_aligned, query_aligned]):
            if aligned[i] != '-':
                if i > 0:
                    sites[j, i] = sites[j, i-1] + 1
                else:
                    sites[j, i] = 1
            else:
                if i > 0:
                    sites[j, i] = sites[j, i-1]
                else:
                    sites[j, i] = 0
    # do not count missing terminal sequences as variants and/or gaps
    gap_start = 0
    gap_end = 0
    for i in range(alignment_length):
        if query_aligned[i] == '-':
            match[i] = 1
        else:
            gap_start = i
            break
    for i in range(alignment_length):
        j = alignment_length - i - 1
        if query_aligned[j] == '-':
            match[j] = 1
        else:
            gap_end = i
            break
    # make variant table
    variants = []
    index = []
    edges = np.where(np.pad(match, [(1,1)], constant_values=1)[1:] != np.pad(match, [(1,1)], constant_values=1)[:-1])[0]
    edges = edges.reshape((edges.shape[0]//2, 2))
    for start, end in zip(edges[:,0], edges[:,1]):
        site_ref = sites[0, start]
        site_con = sites[1, start]
        ref = target_aligned[start:end].replace('-', '')
        alt = query_aligned[start:end].replace('-', '')
        if len(alt) != len(ref) and site_ref > 0:
            # longshot vcf format: left align with one matching base in front of indels
            site_ref -= 1
            site_con -= 1
            ref = target_aligned[start-1:end].replace('-', '')
            alt = query_aligned[start-1:end].replace('-', '')
        variant = "{}{}{}".format(ref, site_ref, alt)
        variants.append( (site_ref, site_con, ref, alt, 'confirmed') )
        index.append(variant)
    var_df = pd.DataFrame(variants, 
                          columns=['site', 'consensus site', 'ref', 'alt', 'decision'], 
                          index=index)
    # make masked regions table
    masked_regions = []
    index = []
    edges = np.where(np.pad(masked, [(1,1)], constant_values=0)[1:] != np.pad(masked, [(1,1)], constant_values=0)[:-1])[0]
    edges = edges.reshape((edges.shape[0]//2, 2))
    for start, end in zip(edges[:,0], edges[:,1]):
        site = sites[0, start]
        masked_ref_bases = len(target_aligned[start:end].replace('-', ''))
        masked_regions.append( (site, site + masked_ref_bases -1, 
                                sites[1, start], sites[1, start] + end - start - 1, 
                                end - start) )
        index.append(site)
    masked_df = pd.DataFrame(masked_regions, 
                             columns=['reference start', 'reference end', 'consensus start', 'consensus end', 'consensus bases'], 
                             index=index)
    return var_df, masked_df, gap_start, gap_end

def get_low_coverage_regions(cov, threshold):
    below = (cov < threshold).astype(np.int16)
    below = np.pad(below, [(1,1)], constant_values=0)
    mask = np.where(below[1:] != below[:-1])[0]
    mask = mask.reshape((mask.shape[0]//2, 2))
    mask[:,0] += 1 # 1-based, inclusive
    df = pd.DataFrame(mask, columns=['start', 'end'])
    df['width'] = df['end'] - df['start'] + 1
    total_bases = np.sum(below)
    return df, total_bases

def dataset_completion_test(sample, artic_runs, nanopore_runs, amplicons):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    # check that (merged) data_dir contains bam files for all primer pools expected to be there based
    # on the nanopore run information
    sample_schemes = nanopore_runs.loc[sample, 'scheme']
    if type(sample_schemes) == str:
        sample_schemes = [sample_schemes]
    else:
        sample_schemes = list(sample_schemes.drop_duplicates())
    for scheme in sample_schemes:
        for pool in list(amplicons[scheme]['pool'].drop_duplicates()):
            bam_fn = os.path.join(artic_dir, "{}.primertrimmed.{}.sorted.bam".format(sample, pool))
            if not os.path.exists(bam_fn):
                logger.warning('Bam file of primer pool {} missing for sample {}'.format(pool, sample))
