import os
import sys
import numpy as np
import pandas as pd
from augur.align import run as augur_align
from Bio import AlignIO, SeqIO

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

def parse_alignment(alignment, samples=None):
    '''parses a Bio.Align.MultipleSeqAlignment object with the reference sequence being at
    index location 0 and returns a pandas dataframe containing variant information
    as well as one containing regions masked with Ns'''
    if not samples:
        samples = [record.id for record in alignment[1:]]
    if type(samples) != list:
        samples = [samples]
    alignment_length = len(alignment[0].seq)
    n_aligned = len(alignment)
    for s, sample in enumerate(samples):
        if len(alignment[s+1].seq) != alignment_length:
            logger.warning('Sample {} seems to be improperly aligned to reference since ' +\
             'its alignment string length does not match the one of the reference'.format(sample))
    assert alignment_length == len(alignment[0].seq)
    sites = np.empty(shape=(n_aligned, alignment_length), dtype=np.int32)
    match = np.empty(shape=(n_aligned, alignment_length), dtype=np.int32)
    masked = np.empty(shape=(n_aligned, alignment_length), dtype=np.int32)
    for i in range(alignment_length):
        for s in range(n_aligned):
            if alignment[s, i].upper() == alignment[0, i].upper():
                if alignment[0, i] == '-':
                    if i > 0:
                        match[s, i] = match[s, i-1]
                        masked[s, i] = masked[s, i-1]
                    else:
                        match[s, i] = 1
                        masked[s, i] = 0
                else:
                    match[s, i] = 1
                    masked[s, i] = 0
            elif alignment[s, i].upper() == 'N' and alignment[0, i].upper() != '-':
                match[s, i] = 1
                masked[s, i] = 1
            else:
                match[s, i] = 0
                masked[s, i] = 0

            if alignment[s, i] != '-':
                if i > 0:
                    sites[s, i] = sites[s, i-1] + 1
                else:
                    sites[s, i] = 1
            else:
                if i > 0:
                    sites[s, i] = sites[s, i-1]
                else:
                    sites[s, i] = 0
    # do not count missing terminal sequences as gaps
    gap_start = {sample:0 for sample in samples}
    gap_end = {sample:0 for sample in samples}
    for s in range(1, n_aligned):
        for i in range(alignment_length):
            if alignment[s, i] == '-':
                match[s, i] = 1
                gap_start[samples[s-1]] += int(alignment[0, i] != '-')
            else:
                break
        for i in range(alignment_length):
            j = alignment_length - i - 1
            if alignment[s, j] == '-':
                match[s, j] = 1
                gap_end[samples[s-1]] += int(alignment[0, j] != '-')
            else:
                break

    variants = []
    multiindex_rows = [[],[]]
    for s, sample in enumerate(samples):
        s += 1
        edges = np.where(np.pad(match[s], [(1,1)], constant_values=1)[1:] != np.pad(match[s], [(1,1)], constant_values=1)[:-1])[0]
        edges = edges.reshape((edges.shape[0]//2, 2))
        for start, end in zip(edges[:,0], edges[:,1]):
            site = sites[0, start]
            ref = str(alignment[0].seq[start:end]).replace('-', '').upper()
            alt = str(alignment[s].seq[start:end]).replace('-', '').upper()
            if len(alt) != len(ref) and site > 0:
                # longshot vcf format: left align with one matching base in front of deletion
                site -= 1
                ref = str(alignment[0].seq[start-1:end]).replace('-', '').upper()
                alt = str(alignment[s].seq[start-1:end]).replace('-', '').upper()
            variant = "{}{}{}".format(ref, site, alt)
            variants.append( (site, ref, alt, 'confirmed') )
            multiindex_rows[0].append(sample)
            multiindex_rows[1].append(variant)
    var_df = pd.DataFrame(variants, 
                          columns=['site', 'ref', 'alt', 'decision'], 
                          index=pd.MultiIndex.from_arrays(multiindex_rows, names=("sample", "variant")))
    masked_regions = []
    multiindex_rows = [[],[]]
    for s, sample in enumerate(samples):
        s += 1
        edges = np.where(np.pad(masked[s], [(1,1)], constant_values=0)[1:] != np.pad(masked[s], [(1,1)], constant_values=0)[:-1])[0]
        edges = edges.reshape((edges.shape[0]//2, 2))
        for start, end in zip(edges[:,0], edges[:,1]):
            site = sites[0, start]
            masked_ref_bases = len(str(alignment[0].seq[start:end]).replace('-', ''))
            masked_regions.append( (site, site + masked_ref_bases -1, 
                                    sites[s, start], sites[s, start] + masked_ref_bases -1, 
                                    masked_ref_bases) )
            multiindex_rows[0].append(sample)
            multiindex_rows[1].append(site)
    masked_df = pd.DataFrame(masked_regions, 
                             columns=['reference start', 'reference end', 'consensus start', 'consensus end', 'bases'], 
                             index=pd.MultiIndex.from_arrays(multiindex_rows, names=("sample", "site")))
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
