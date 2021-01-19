import os
import sys
import re
from glob import glob
import numpy as np
import pandas as pd
from jinja2 import Environment, PackageLoader, select_autoescape
from .verify import dataset_completion_test, get_masked_bases, get_low_coverage_regions
from .data_parser import *
from .plotting import create_summary_plot

import logging
logger = logging.getLogger()

pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option('display.width', 130)

def sample_report(sample, template, sample_results_dir, sample_schemes, cov_primertrimmed, 
                  cov_illumina, cov_sanger, cov_pools, snv_info, reference, reference_genes, 
                  amplicons, sample_nanopore_runs, threshold_limit, threshold_low, 
                  filtered_snvs_only=True):
    img_dir = os.path.join(sample_results_dir, 'img')
    report_fn = os.path.join(sample_results_dir, "{}.cov2seq-report.html".format(sample))
    summary_plot_fn = os.path.join(img_dir, "{}.report.svg".format(sample))
    final_consensus_fn = os.path.join(sample_results_dir, "{0}.final.fasta".format(sample))
    cov_limit_regions, cov_limit_bases = get_low_coverage_regions(cov_primertrimmed, threshold_limit)
    cov_low_regions, cov_low_bases = get_low_coverage_regions(cov_primertrimmed, threshold_low)
    final = os.path.exists(final_consensus_fn)

    header = "Report for sample {}".format(sample)
    logger.debug('final: {}'.format(final))
    if not final:
        header += " (preliminary)"
    if filtered_snvs_only:
        snv_info = snv_info.loc[pd.notnull(snv_info[('ARTIC', 'snv_filter')]) |\
                                pd.notnull(snv_info[('final', 'decision')])]
    snv_info_ = snv_info.reset_index().reset_index().set_index(['index', 'level_0'])
    snv_info_.index.set_names(['variant', ''], inplace=True)
    snv_table = snv_info_.to_html(float_format=lambda f: "{:.1f}".format(f),
                                  na_rep="", justify="left")
    nanopore_runs_table = sample_nanopore_runs.reset_index(drop=True).to_html()
    render_dict = {"sample" : sample,
                   "header" : header,
                   "snv_table" : snv_table,
                   "summary_plot_fn" : summary_plot_fn,
                   "threshold_limit" : threshold_limit,
                   "threshold_low" : threshold_low,
                   "nanopore_runs_table" : nanopore_runs_table}
    # create output directory for sample if it does not exist already
    dirs = [sample_results_dir, img_dir]
    for d in dirs:
        if not os.path.exists(d):
            os.makedirs(d)
            cmd = 'chmod g+w {}'.format(d)
            os.system(cmd)
        if not os.access(d, os.W_OK):
            logger.error('Missing write permission for directory {}'.format(d))
            exit(1)

    create_summary_plot(sample_schemes, cov_primertrimmed, cov_illumina, cov_sanger, cov_pools, 
                        cov_limit_regions, cov_low_regions, snv_info, reference, reference_genes,
                        amplicons, final, savepath=summary_plot_fn)

    with open(report_fn, 'w') as f:
        print(template.render(**render_dict), file=f)

def selected_samples_from_patterns(nanopore_dir, patterns):
    sample_paths = []
    for pattern in patterns:
        sample_paths.extend(glob(os.path.join(nanopore_dir, pattern)))
    selected_samples = [os.path.basename(p)  for p in set(sample_paths)]
    selected_samples.sort()
    return selected_samples

def create_sample_reports(args, pkg_dir):
    selected_samples = selected_samples_from_patterns(args.nanopore_dir, args.samples)
    logger.info("Creating reports for the following {} samples:\n{}".format(
        len(selected_samples), ", ".join(selected_samples)))
    reference, reference_genes = load_reference(args.reference_fn, args.reference_annotation_fn)
    clades_df, subclades_df = load_clades_info(args.nextstrain_ncov)
    artic_runs, nanopore_runs = load_nanopore_info(selected_samples, args.nanopore_dir)
    primer_schemes = list(nanopore_runs.scheme.drop_duplicates())
    primers, amplicons = load_primer_schemes(args.primer_schemes_dir, primer_schemes)

    logger.info(pkg_dir)
    jinja_env = Environment(
        loader=PackageLoader('cov2seq', 'templates'),
        autoescape=False
    )
    template = jinja_env.get_template('report.html')

    for sample in selected_samples:
        logger.info('Creating sample report for sample {}'.format(sample))
        dataset_completion_test(sample, artic_runs, nanopore_runs, amplicons)

        sample_results_dir = os.path.join(args.results_dir, sample)
        sample_schemes = nanopore_runs.loc[sample, 'scheme']
        if type(sample_schemes) == str:
            sample_schemes = [sample_schemes]
        else:
            sample_schemes = list(sample_schemes.drop_duplicates())
        sample_nanopore_runs = nanopore_runs.loc[sample]
        if type(sample_nanopore_runs) == pd.core.series.Series:
            sample_nanopore_runs = sample_nanopore_runs.to_frame().T
        cov_primertrimmed = get_nanopore_coverage_primertrimmed(sample, artic_runs)
        cov_illumina, mapped_illumina = get_illumina_coverage_and_mappings(sample, args.illumina_dir, reference)
        cov_sanger = approximate_sanger_coverage(sample, args.sanger_dir, reference, amplicons, primers)
        cov_pools = get_nanopore_pool_coverage(sample, artic_runs, nanopore_runs, amplicons, reference)
        snv_info = load_snv_info(sample, artic_runs, args.results_dir, args.reference_fn, args.reference_annotation_fn, clades_df)

        sample_report(sample, template, sample_results_dir, sample_schemes, cov_primertrimmed, 
                      cov_illumina, cov_sanger, cov_pools, snv_info, reference, reference_genes, 
                      amplicons, sample_nanopore_runs, args.threshold_limit, args.threshold_low)
