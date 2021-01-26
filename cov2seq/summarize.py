import os
from .data_parser import *

def create_resequencing_scheme(args):
    selected_samples = selected_samples_from_patterns(args.nanopore_dir, args.samples)
    if len(selected_samples) == 0:
        logger.error('No samples found in {} matching the given pattern(s)'.format(args.nanopore_dir))
        exit(1)
    logger.info("Creating a re-sequencing scheme for the following {} samples:\n{}".format(
        len(selected_samples), ", ".join(selected_samples)))

    artic_runs, nanopore_runs = load_nanopore_info(selected_samples, args.nanopore_dir,
                                                   sorted_reads=False, trimmed_reads=False, primertrimmed_reads=False)
    primer_schemes = list(nanopore_runs.scheme.drop_duplicates())
    primers, amplicons = load_primer_schemes(args.primer_schemes_dir, primer_schemes)