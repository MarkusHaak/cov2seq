import os
from .data_parser import *
from .verify import get_low_coverage_regions
from .plotting import plot_resequencing_scheme

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatch

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def get_amplicons_by_regions(regions, scheme, amplicons):
    affected_amplicons = {}
    overlap_only = []
    for i,(start,end,bases) in regions.iterrows():
        j = 0
        for i,amplicon in amplicons[scheme].iterrows():
            if j > 0:
                prv_amplicon = amplicons[scheme].iloc[j-1]
            else:
                prv_amplicon = {'end' : amplicon['start']}
            if j < amplicons[scheme].shape[0] - 1:
                nxt_amplicon = amplicons[scheme].iloc[j+1]
            else:
                nxt_amplicon = {'start' : amplicon['end']}
                            
            if amplicon['start'] <= start < nxt_amplicon['start']:
                intersection = min(end, amplicon['end']) - start + 1
                if amplicon['name'] not in affected_amplicons:
                    affected_amplicons[amplicon['name']] = [1, intersection] # #regions, #bases
                else:
                    affected_amplicons[amplicon['name']][0] += 1
                    affected_amplicons[amplicon['name']][1] += intersection
            elif prv_amplicon['end'] < end <= amplicon['end']:
                intersection =  end - max(start, amplicon['start']) + 1
                if amplicon['name'] not in affected_amplicons:
                    affected_amplicons[amplicon['name']] = [1, intersection] # #regions, #bases
                else:
                    affected_amplicons[amplicon['name']][0] += 1
                    affected_amplicons[amplicon['name']][1] += intersection
            elif start < amplicon['start'] < amplicon['end'] < end:
                intersection = amplicon['end'] - amplicon['start'] + 1
                if amplicon['name'] not in affected_amplicons:
                    affected_amplicons[amplicon['name']] = [1, intersection] # #regions, #bases
                else:
                    affected_amplicons[amplicon['name']][0] += 1
                    affected_amplicons[amplicon['name']][1] += intersection
            elif nxt_amplicon['start'] <= start <= end <= amplicon['end']:
                overlap_only.append( (start,end) )
            j += 1
    for start, end in overlap_only:
        j = 0
        for i,amplicon in amplicons[scheme].iterrows():
            if j < amplicons[scheme].shape[0] - 1:
                nxt_amplicon = amplicons[scheme].iloc[j+1]
            else:
                nxt_amplicon = {'start' : amplicon['end'], 'name' : amplicon['name']}
            if nxt_amplicon['start'] <= start <= end <= amplicon['end']:
                if nxt_amplicon['name'] in affected_amplicons:
                    affected_amplicons[nxt_amplicon['name']][0] += 1
                    affected_amplicons[nxt_amplicon['name']][1] += intersection
                else:
                    affected_amplicons[amplicon['name']][0] += 1
                    affected_amplicons[amplicon['name']][1] += intersection
    return affected_amplicons

def create_resequencing_scheme(args):
    selected_samples = selected_samples_from_patterns(args.nanopore_dir, args.samples)
    if len(selected_samples) == 0:
        logger.error('No samples found in {} matching the given pattern(s)'.format(args.nanopore_dir))
        exit(1)
    logger.info("Creating a re-sequencing scheme for the following {} samples:\n{}".format(
        len(selected_samples), ", ".join(selected_samples)))
    #reference_fn = os.path.join(args.nextstrain_ncov_dir, "defaults", "reference_seq.gb")
    #reference, reference_genes = load_reference(reference_fn)
    artic_runs, nanopore_runs = load_nanopore_info(selected_samples, args.nanopore_dir,
                                                   sorted_reads=False, trimmed_reads=False, primertrimmed_reads=False)
    primer_schemes = list(nanopore_runs.scheme.drop_duplicates())
    primers, amplicons = load_primer_schemes(args.primer_schemes_dir, primer_schemes)

    for scheme in primer_schemes:
        mat = np.zeros(shape=(len(selected_samples),len(amplicons[scheme])), dtype=np.float32)
        for i,sample in enumerate(selected_samples):
            cov_primertrimmed = get_nanopore_coverage_primertrimmed(sample, artic_runs)
            
            limit_regions,_ = low_regions,_ = get_low_coverage_regions(cov_primertrimmed, args.threshold_limit)
            limit_affected_amplicons = get_amplicons_by_regions(limit_regions, scheme, amplicons)
            low_regions,_ = get_low_coverage_regions(cov_primertrimmed, args.threshold_low)
            #low_affected_amplicons = get_amplicons_by_regions(low_regions, scheme, amplicons)
            low_regions_filtered = low_regions.loc[low_regions.width > args.tolerated_consecutive_bases]
            low_affected_amplicons_filtered = get_amplicons_by_regions(low_regions_filtered, scheme, amplicons)
            j = 0
            for amplicon,(_) in amplicons[scheme].iterrows():
                if amplicon in limit_affected_amplicons:
                    mat[i,j] = 2
                elif amplicon in low_affected_amplicons_filtered:
                    mat[i,j] = 1
                else:
                    mat[i,j] = 0
                j += 1
        plot_resequencing_scheme(selected_samples, amplicons[scheme], mat, 
                                 figwidth=args.figwidth, max_cols_per_pool=args.max_cols_per_pool, 
                                 fontsize=args.fontsize)

