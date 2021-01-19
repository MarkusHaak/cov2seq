import os
import sys
import numpy as np
import pandas as pd
from augur.align import run as augur_align
from Bio import AlignIO, SeqIO

import logging
logger = logging.getLogger()

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
    return alignment

def variants_from_alignment(alignment, samples=None):
    '''translates a Bio.Align.MultipleSeqAlignment object with the reference sequence being at
    index location 0 to a pandas dataframe containing variant information'''
    if not samples:
        samples = [record.id for record in alignment[1:]]
    if type(samples) != list:
        samples = [samples]
    alignment_length = len(alignment[0].seq)
    n_aligned = len(alignment)
    for s, sample in enumerate(samples):
        if len(alignment[s+1].seq) != alignment_length:
            logger.warning('Sample {} seems to be improperly aligned to reference since its alignment string length does not match the one of the reference'.format(sample))
    assert alignment_length == len(alignment[0].seq)
    locs = np.empty(shape=(len(alignment), alignment_length), dtype=np.int32)
    match = np.empty(shape=(len(alignment), alignment_length), dtype=np.int32)
    for i in range(alignment_length):
        for s in range(n_aligned):
            if int((alignment[s, i].upper() == alignment[0, i].upper() or alignment[s, i].upper() == 'N') and alignment[0, i] != '-'):
                match[s, i] = 1
            elif alignment[s, i] == '-' and alignment[0, i] != '-':
                if i > 0:
                    match[s, i] = match[s, i-1]
                else:
                    match[s, i] = 1
            else:
                match[s, i] = 0
            if alignment[s, i] != '-':
                if i > 0:
                    locs[s, i] = locs[s, i-1] + 1
                else:
                    locs[s, i] = 1
    variants = []
    multiindex_rows = [[],[]]
    for s, sample in enumerate(samples):
        s += 1
        edges = np.where(np.pad(match[s,1:], [(1,1)], constant_values=1) != np.pad(match[s,:-1], [(1,1)], constant_values=1))[0]
        edges = edges.reshape((edges.shape[0]//2, 2))
        for start, end in zip(edges[:,0], edges[:,1]):
            site = locs[0, start]
            ref = str(alignment[0].seq[start:end]).replace('-', '')
            alt = str(alignment[s].seq[start:end]).replace('-', '')
            variant = "{}{}{}".format(ref, site, alt)
            variants.append( (site, ref, alt) )
            multiindex_rows[0].append(sample)
            multiindex_rows[1].append(variant)
    df = pd.DataFrame(variants, columns=['site', 'REF', 'ALT'], index=pd.MultiIndex.from_arrays(multiindex_rows, names=("sample", "variant")))
    return df

def get_masked_bases(fasta_fn):
    record = next(SeqIO.parse(fasta_fn, 'fasta'))
    n_sites = np.array([letter.upper() == 'N' for letter in record.seq], dtype=np.int16)
    n_sites = np.pad(n_sites, [(1,1)], constant_values=0)
    mask = np.where(n_sites[1:] != n_sites[:-1])[0]
    mask = mask.reshape((mask.shape[0]//2, 2))
    mask[:,0] += 1 # 1-based, inclusive
    df = pd.DataFrame(mask, columns=['start', 'end'])
    df['width'] = df['end'] - df['start'] + 1
    total_bases = np.sum(n_sites)
    return df, total_bases

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
