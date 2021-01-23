import os
import sys
import re
from glob import glob
import numpy as np
import pandas as pd
from Bio import SeqIO
from jinja2 import Environment, PackageLoader, select_autoescape
import pkg_resources
from .verify import dataset_completion_test, get_low_coverage_regions
from .data_parser import *
from .plotting import create_summary_plot

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

pd.set_option("display.max_rows", None, "display.max_columns", None)
pd.set_option('display.width', 130)

def txt_color(s, color):
    return f'<div style="color:{color}">{s}</div>'

def export_vcf(sample, snv_info, sample_results_dir, reference):
    vcf_fn = os.path.join(sample_results_dir, '{}.final.vcf'.format(sample))
    logger.info('Exporting vcf file {}'.format(vcf_fn))
    with open(vcf_fn, 'w') as f:
        print('##fileformat=VCFv4.2', file=f)
        print('##source=cov2seq v{}'.format(pkg_resources.get_distribution('cov2seq').version), file=f)
        print('##reference={}'.format(reference.id), file=f)
        print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', file=f)
        sel = snv_info[('final', 'decision')].isin(['confirmed', 'introduced']) & \
              ~snv_info.index.duplicated(keep='first')
        for i,row in snv_info.loc[sel].iterrows():
            CHROM = reference.id
            POS = str(row[('medaka variant', 'site')])
            ID = '.'
            REF = row[('medaka variant', 'ref')]
            ALT = row[('medaka variant', 'alt')]
            QUAL = '.'
            FILTER = 'PASS'
            INFO = '.'
            print("\t".join([CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO]), file=f)

    if os.path.exists(vcf_fn):
        cmd = 'chmod g+w {}'.format(vcf_fn)
        os.system(cmd)

def sample_report(sample, template, sample_results_dir, sample_schemes, cov_primertrimmed, 
                  cov_illumina, cov_sanger, cov_pools, snv_info, reference, reference_genes, 
                  amplicons, sample_nanopore_runs, sample_artic_stats, clade_assignment, parent_clade,
                  masked_regions, threshold_limit, threshold_low, software_versions, 
                  filtered_snvs_only=True):
    img_dir = os.path.join(sample_results_dir, 'img')
    report_fn = os.path.join(sample_results_dir, "{}.report.html".format(sample))
    summary_plot_fn = os.path.join(img_dir, "{}.summary.svg".format(sample))
    final_consensus_fn = os.path.join(sample_results_dir, "{0}.final.fasta".format(sample))
    cov_limit_regions, cov_limit_bases = get_low_coverage_regions(cov_primertrimmed, threshold_limit)
    cov_low_regions, cov_low_bases = get_low_coverage_regions(cov_primertrimmed, threshold_low)
    final = os.path.exists(final_consensus_fn)
    decision_colors = {'confirmed' : 'green',
                       'masked': 'gray',
                       'introduced': 'magenta',
                       'partially masked': 'yellow'}
    impact_colors = {'LOW' : 'black',
                     'MODIFIER': 'blue',
                     'MODERATE': 'orange',
                     'HIGH': 'red'}
    coverage_colors = {i : 'red' if i < threshold_limit else 'orange' for i in range(threshold_low)}

    header = "Report for sample {}".format(sample)
    logger.debug('final: {}'.format(final))
    if not final:
        header += " (preliminary)"
    snv_info_ = snv_info
    if filtered_snvs_only:
        snv_info_ = snv_info.loc[pd.notnull(snv_info[('ARTIC', 'snv_filter')]) |\
                                 pd.notnull(snv_info[('final', 'decision')])]
    snv_info_ = snv_info_.reset_index().reset_index().set_index(['index', 'level_0'])
    snv_info_.index.set_names(['variant', ''], inplace=True)
    formatters = {
        ('longshot', 'qual') : lambda x: txt_color("{:.1f}".format(x), "orange") if x<500. else "{:.1f}".format(x),
        ('longshot', 'cov') : lambda x: txt_color(int(x), coverage_colors.get(int(x), 'black')),
        ('longshot', '#ref') : lambda x: "{:.0f}".format(x),
        ('longshot', '#alt') : lambda x: "{:.0f}".format(x),
        ('longshot', '#amb') : lambda x: "{:.0f}".format(x),
        ('longshot', 'strand bias') : lambda x: txt_color(x, 'orange') if x==True else f"{x}",
        ('ARTIC', 'snv_filter') : lambda x: txt_color(x, 'red') if x==False else txt_color(x, 'green'),
        ('final', 'decision') : lambda x: txt_color(x, decision_colors.get(x, 'red')),
        ('snpEff', 'impact') : lambda x: txt_color(x, impact_colors.get(x, 'black'))
    }
    snv_table = snv_info_.to_html(float_format=lambda f: "{:.1f}".format(f),
                                  na_rep="", justify="left", classes=['table-hover'],
                                  formatters=formatters, escape=False)
    nanopore_runs_table = sample_nanopore_runs.reset_index(drop=True).to_html(classes=['table-hover'])
    artic_stats_table = sample_artic_stats.to_html(classes=['table-hover'])
    software_versions_table = software_versions.to_html()
    if final:
        masked_regions_table = masked_regions.to_html(classes=['table-hover'])
        masked_bases = np.sum(masked_regions['bases'])
        consensus_length = len(next(SeqIO.parse(final_consensus_fn, "fasta")).seq)
    else:
        masked_bases, masked_regions_table, consensus_length = None, None, None
    render_dict = {"sample" : sample,
                   "final": final,
                   "header" : header,
                   "snv_table" : snv_table,
                   "summary_plot_fn" : summary_plot_fn,
                   "threshold_limit" : threshold_limit,
                   "threshold_low" : threshold_low,
                   "nanopore_runs_table" : nanopore_runs_table,
                   "artic_stats_table" : artic_stats_table,
                   "medaka_variant_SNVs": len(snv_info),
                   "unique_medaka_variant_SNVs": len(snv_info.index.unique()),
                   "longshot_SNVs": len(snv_info_.index.droplevel(1).unique()),
                   "clade_assignment": clade_assignment, 
                   "parent_clade": parent_clade,
                   "software_versions_table": software_versions_table,
                   "masked_bases": masked_bases,
                   "masked_regions_table": masked_regions_table,
                   "consensus_length": consensus_length}
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

    for fn in [report_fn, summary_plot_fn]:
        if os.path.exists(fn):
            cmd = 'chmod g+w {}'.format(fn)
            os.system(cmd)

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
    reference_fn = os.path.join(args.nextstrain_ncov_dir, "defaults", "reference_seq.gb")
    reference_fasta_fn = os.path.join(args.nextstrain_ncov_dir, "defaults", "reference_seq.fasta")
    reference, reference_genes = load_reference(reference_fn)
    clades_df, subclades_df = load_clades_info(args.nextstrain_ncov_dir)
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
        final_consensus_fn = os.path.join(sample_results_dir, "{0}.final.fasta".format(sample))
        sample_schemes = nanopore_runs.loc[sample, 'scheme']
        if type(sample_schemes) == str:
            sample_schemes = [sample_schemes]
        else:
            sample_schemes = list(sample_schemes.drop_duplicates())
        sample_nanopore_runs = nanopore_runs.loc[sample]
        if type(sample_nanopore_runs) == pd.core.series.Series:
            sample_nanopore_runs = sample_nanopore_runs.to_frame().T
        sample_artic_stats = artic_runs.loc[sample].to_frame().T
        cov_primertrimmed = get_nanopore_coverage_primertrimmed(sample, artic_runs)
        cov_illumina, mapped_illumina = get_illumina_coverage_and_mappings(sample, args.illumina_dir, reference)
        cov_sanger = approximate_sanger_coverage(sample, args.sanger_dir, reference, amplicons, primers)
        cov_pools = get_nanopore_pool_coverage(sample, artic_runs, nanopore_runs, amplicons, reference)
        snv_info, masked_regions = load_snv_info(sample, artic_runs, args.results_dir, 
                                                 reference_fasta_fn, args.snpeff_dir, clades_df, 
                                                 subclades_df)
        _,clade_assignment, parent_clade = assign_clade(sample, artic_runs, args.results_dir, 
                                                      args.nextstrain_ncov_dir, repeat_assignment=True)
        software_versions = get_software_versions(sample, artic_runs, args.results_dir)

        sample_report(sample, template, sample_results_dir, sample_schemes, cov_primertrimmed, 
                      cov_illumina, cov_sanger, cov_pools, snv_info, reference, reference_genes, 
                      amplicons, sample_nanopore_runs, sample_artic_stats, clade_assignment, parent_clade,
                      masked_regions, args.threshold_limit, args.threshold_low, software_versions)

        if args.export_vcf and os.path.exists(final_consensus_fn):
            export_vcf(sample, snv_info, sample_results_dir, reference)
