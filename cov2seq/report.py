import os
import sys
from glob import glob
import numpy as np
import pandas as pd
from .verify import dataset_completion_test, get_masked_bases
from .data_parser import *

import logging
logger = logging.getLogger()

def get_selected_samples(nanopore_dir, patterns):
    sample_paths = []
    for pattern in patterns:
        sample_paths.extend(glob(os.path.join(nanopore_dir, pattern)))
    selected_samples = [os.path.basename(p)  for p in set(sample_paths)]
    selected_samples.sort()
    return selected_samples

def create_sample_reports(args):
    selected_samples = get_selected_samples(args.nanopore_dir, args.samples)
    logger.info("Creating reports for the following {} samples:\n{}".format(
        len(selected_samples), ", ".join(selected_samples)))
    reference, reference_genes = load_reference(args.reference_fn, args.reference_annotation_fn)
    clades_df = load_clades_info(args.nextstrain_ncov)
    artic_runs, nanopore_runs = load_nanopore_info(selected_samples, args.nanopore_dir)
    primer_schemes = list(nanopore_runs.scheme.drop_duplicates())
    primers, amplicons = load_primer_schemes(args.primer_schemes_dir, primer_schemes)
    for sample in selected_samples:
        logger.info(logger.info('Creating sample report for sample {}'.format(sample)))
        dataset_completion_test(sample, artic_runs, nanopore_runs, amplicons)
        
