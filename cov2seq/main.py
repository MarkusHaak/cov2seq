import os
import sys
import configparser
import pathlib
import pkg_resources
import argparse
from .analyze import run_artic_medaka
from .report import create_sample_reports
from .summarize import create_resequencing_scheme

import logging
logger = logging.getLogger()
shandler = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s: %(message)s')
shandler.setFormatter(formatter)
logger.addHandler(shandler)
logger.setLevel(logging.INFO)

class ArgHelpFormatter(argparse.HelpFormatter):
    '''
    Formatter adding default values to help texts.
    '''
    def __init__(self, prog):
        super().__init__(prog)

    def _get_help_string(self, action):
        text = action.help
        if action.default is not None and \
           action.default != argparse.SUPPRESS and \
           '(default' not in text.lower():
           text += ' ({})'.format(action.default)
        return text

def get_package_info():
    pkg_dir = os.path.dirname(pathlib.Path(__file__).parent.absolute())
    pkg_res = pkg_resources.get_distribution('cov2seq')
    pkg_version = pkg_res.version
    pkg_descr = list(pkg_res._get_metadata(pkg_res.PKG_INFO))[-2]
    return pkg_dir, pkg_version, pkg_descr

def read_configuration(argv, pkg_dir, defaults=None, sections=None):
    def update_defaults(conf_file):
        config = configparser.SafeConfigParser()
        config.read([conf_file])
        defaults.update(dict(config.items('DEFAULT')))
        for section in sections:
            if config.has_section(section):
                defaults.update(dict(config.items(section)))

    if sections is None:
        sections = []
    if defaults is None:
        defaults = {}
    default_cfg = os.path.join(pkg_dir, "cov2seq", "cov2seq.cfg")
    user_cfg = os.path.join(os.path.expanduser('~'), "cov2seq.cfg")
    local_cfg = "cov2seq.cfg"

    conf_parser = argparse.ArgumentParser(
        description=__doc__,        
        formatter_class=argparse.RawDescriptionHelpFormatter,
        add_help=False)
    conf_parser.add_argument("-c", "--conf_file",
        help='''Specify an additional config file location with highest priority. 
             (defaults, increasing priority: {})'''.format([default_cfg, user_cfg, local_cfg]), 
        metavar="FILE")
    args, remaining_argv = conf_parser.parse_known_args()

    if os.path.exists(default_cfg):
        update_defaults(default_cfg)
    else:
        logger.warning("Package-default configuration file not found.")
    if os.path.exists(user_cfg):
        update_defaults(user_cfg)
    else:
        logger.debug("No user configuration file found at {}.".format(user_cfg))
    if os.path.exists(local_cfg):
        update_defaults(local_cfg)
    else:
        logger.debug("No local configuration file {} found.".format(local_cfg))
    if args.conf_file:
        update_defaults(args.conf_file)

    return conf_parser, remaining_argv, defaults

def add_main_group_to_parser(parser):
    main_group = parser.add_argument_group('Main Options')
    main_group.add_argument('--nanopore_dir',
        help='''base directory containing sample-named sub-directories of nanopore datasets analyzed
             with the artic-ncov2019 pipeline. The directory tree must be 
             <nanopore_dir>/<sample_name>/artic-medaka/<scheme_name>/<scheme_version>/<sample_name>.*''')
    main_group.add_argument('--illumina_dir',
        help='''base directory containing sample-named sub-directories of illumina datasets. 
             The directory tree must be <illumina_dir>/<sample_name>/*.bam''')
    main_group.add_argument('--sanger_dir',
        help='''base directory containing sample-named sub-directories of sanger data. 
             The directory tree must be <sanger_dir>/<sample_name>/*-<primer_id>.ab1''')
    main_group.add_argument('--primer_schemes_dir',
        help='''Direcotry containing the primer scheme files in a directory tree as follows:
             <primer_schemes_dir>/<scheme_name>/<scheme_version>/<scheme_name>.*''')
    main_group.add_argument('--results_dir',
        help='''Directory for saving sample results.''')
    main_group.add_argument('--nextstrain_ncov_dir',
        help='''Path to the directory containing the Nextstrain ncov git repository. This should
             be updated regularly to get the latest nextstrain clade information.''')
    main_group.add_argument('--snpeff_dir',
        help='''Path to the directory containing the snpEff binaries and the snpEff data directory.''')
    main_group.add_argument('--ct_values_fn',
        help='''Csv file containing sample names in the first column and their ct values in the second 
             column (semicolon column separator, header in first line).''')
    main_group.add_argument('--threshold_limit',
        help='''Minimal acceptable nanopore coverage. Regions with a coverage below this 
                this threshold are highlighted red.''',
        type=int)
    main_group.add_argument('--threshold_low',
        help='''Minimal acceptable nanopore coverage to assume that variant calling based on nanopore 
                data alone works sufficiently well. Regions with a coverage below this 
                this threshold are highlighted orange.''',
        type=int)
    return parser

def check_arguments(args):
    for arg, val in [('nanopore_dir', args.nanopore_dir), ('primer_schemes_dir', args.primer_schemes_dir), 
                     ('results_dir', args.results_dir), ('nextstrain_ncov_dir', args.nextstrain_ncov_dir), 
                     ('snpeff_dir', args.snpeff_dir), ('illumina_dir', args.illumina_dir),
                     ('sanger_dir', args.sanger_dir), ('ct_values_fn', args.ct_values_fn)]:
        if val is None:
            logger.error('No value given for required argument: {}'.format(arg))
            exit(1)
    args.nanopore_dir = os.path.expanduser(args.nanopore_dir)
    args.primer_schemes_dir = os.path.expanduser(args.primer_schemes_dir)
    args.results_dir = os.path.expanduser(args.results_dir)
    args.nextstrain_ncov_dir = os.path.expanduser(args.nextstrain_ncov_dir)
    args.snpeff_dir = os.path.expanduser(args.snpeff_dir)
    args.illumina_dir = os.path.expanduser(args.illumina_dir)
    args.sanger_dir = os.path.expanduser(args.sanger_dir)
    args.ct_values_fn = os.path.expanduser(args.ct_values_fn)
    for arg in [args.nanopore_dir,
                args.results_dir,
                args.primer_schemes_dir,
                args.nextstrain_ncov_dir,
                args.snpeff_dir]:
        if not os.path.isdir(arg):
            logger.error('Directory {} does not exist but is required.'.format(arg))
            exit(1)
        if not os.access(arg, os.W_OK) and arg in [args.nanopore_dir, args.results_dir]:
            logger.error('Write permissions required for directory {}'.format(arg))
            exit(1)
    for arg in [args.illumina_dir,
                args.sanger_dir]:
        if not os.path.isdir(arg):
            logger.warning('Directory {} does not exist.'.format(arg))
    return args

def add_help_group_to_parser(parser):
    pkg_dir, pkg_version, pkg_descr = get_package_info()
    help_group = parser.add_argument_group('Help')
    help_group.add_argument('--version', 
        action='version',
        default=argparse.SUPPRESS,
        version="cov2seq {}".format(pkg_version))
    help_group.add_argument('-h', '--help', 
        action='help', 
        default=argparse.SUPPRESS,
        help='Show this help message and exit.')
    return parser

def init_parser(argv, defaults, script_descr="", sections=[]):
    pkg_dir, pkg_version, pkg_descr = get_package_info()
    conf_parser, remaining_argv, defaults = read_configuration(argv, pkg_dir, defaults, sections)
    descr = pkg_descr
    if script_descr:
        descr += "\n" + script_descr
    parser = argparse.ArgumentParser(
        description=descr,
        parents=[conf_parser],
        formatter_class=ArgHelpFormatter, 
        add_help=False)
    parser.set_defaults(**defaults)
    return add_main_group_to_parser(parser), remaining_argv

def analyze(argv=None):
    if argv is None:
        argv = sys.argv
    defaults = {'restrict': [],
                'exclude': []} # configuration file independent default values; lowest priority
    script_descr = 'This script automates the execution of ARTIC guppyplex and the ARTIC medaka minion pipeline ' +\
                   'for a set of sequencing runs.'
    parser, remaining_argv = init_parser(argv, defaults, script_descr=script_descr, sections=['ANALYZE'])

    analyze_group = parser.add_argument_group('Analyze Option Group')
    analyze_group.add_argument('-i', '--input-directories',
        help='''Paths to nanopore sequencing run directories of SARS-CoV2 samples. Each directory
             must contain contain a fastq_pass directory with demultiplexed fastq files sorted into individual
             barcode subdirectories and the run_configuration.json and primers.json configuration files 
             used for ARTIC rampart.''',
        nargs='+',
        required=True)
    analyze_group.add_argument('-r', '--restrict',
        help='Restrict analysis to the specified sample(s).',
        metavar='sample_id', nargs='+')
    analyze_group.add_argument('-e', '--exclude', 
        help='''Exclude the specified sample(s) from the analysis. Samples with "ntc" in their sample IDs 
             are excluded automatically.''',
        metavar='sample_id', nargs='+')
    analyze_group.add_argument('--dry-run',
        help='only print what would be done without actually executing any commands.', 
        action='store_true')
    analyze_group.add_argument('--overwrite', 
        help='''Overwrite all previous analysis results. Required if a sample directory under <nanopore_dir>,
             for example if an analysis is repeated or a sample was re-sequenced. In the latter case, ALL run 
             directories must be given as input directories or else they will not contribute to the sample's 
             analysis results.''',
        action='store_true')
    analyze_group.add_argument('--overview',
        help='''Parse and print all information present in the configuration files of given input directories
             without executing an analysis.''',
        action='store_true')
    analyze_group.add_argument('--default_scheme',
        halp='''Default primer scheme that is assumed if no primers.json configuration file exists in an input directory.''')
    analyze_group.add_argument('--normalize',
        help='''Normalize coverage to approximately this amount of reads per stand.''',
        type=int)
    analyze_group.add_argument('--min_len',
        help='''Minimum length of reads to be considered in the analysis pipeline. At the moment, this parameter
             cannot be set for every Amplicon set individually.''',
        type=int)
    analyze_group.add_argument('--max_len',
        help='''Maximum length of reads to be considered in the analysis pipeline. At the moment, this parameter
             cannot be set for every Amplicon set individually.''',
        type=int)
    analyze_group.add_argument('--min_quality',
        help='Minimum quality of reads to be considered in the analysis pipeline.',
        type=int)
    analyze_group.add_argument('-t', '--threads', type=int, help='number of threads', default=6)

    parser = add_help_group_to_parser(parser)
    args = parser.parse_args(remaining_argv)
    args = check_arguments(args)
    run_artic_medaka(args)

def summarize(argv=None):
    if argv is None:
        argv = sys.argv
    defaults = {} # configuration file independent default values; lowest priority
    script_descr = "Script for creating figures and tables summarizing and comparing data of multiple samples."
    parser, remaining_argv = init_parser(argv, defaults, script_descr=script_descr, sections=['SUMMARIZE'])

    summarize_group = parser.add_argument_group('Summarize Option Group')
    summarize_group.add_argument('-t', '--type',
        help='Choose which type of summary shall be created.',
        required=True,
        choices=['resequencing-scheme'])
    summarize_group.add_argument('-s', '--samples',
        help='''Glob patterns matching sample names present in <nanopore_dir> directory.
             See https://docs.python.org/3/library/glob.html for a description of valid patterns.''',
        nargs='+',
        required=True)
    summarize_group.add_argument('--tolerated_consecutive_bases',
        help='''Number of tolerated consecutive bases with coverage below <threshold_low> 
             but above <threshold_limit> so as to not create the need for re-sequencing.''',
        type=int)
    summarize_group.add_argument('--scheme',
        help='''Enforces this scheme to be used for re-sequencing. If none is given, 
             an individual re-sequencing schedule is created for every schemes that was so far used for
             sequencing any of the selected samples.''')

    cosmetic_group = parser.add_argument_group('Cosmetic Options Group')
    cosmetic_group.add_argument('--figwidth',
        help='''Width of the figure in inches.''',
        type=int)
    cosmetic_group.add_argument('--max_cols_per_pool',
        help='''Maximum columns of per primer pool for horizontal side-by-side presentation.''',
        type=int)
    cosmetic_group.add_argument('--fontsize',
        help='''Font size.''',
        type=int)

    parser = add_help_group_to_parser(parser)
    args = parser.parse_args(remaining_argv)
    args = check_arguments(args)
    if args.type == 'resequencing-scheme':
        create_resequencing_scheme(args)

def report(argv=None):
    if argv is None:
        argv = sys.argv
    defaults = {} # configuration file independent default values; lowest priority
    script_descr = "This script creates comprehensive reports for one or several samples."
    parser, remaining_argv = init_parser(argv, defaults, script_descr=script_descr, sections=['REPORT'])

    report_group = parser.add_argument_group('Report Option Group')
    report_group.add_argument('-s', '--samples',
        help='Relative glob patterns matching sample names present in <nanopore_dir> directory.',
        nargs='+',
        required=True)
    report_group.add_argument('--repeat_assignment',
        help='Repeat the clade assignment even if a clade assignment was already once performed.',
        action='store_true')
    report_group.add_argument('--alignment_tool',
        help='''Alignment tool to use for pairwise alignment of consensus sequences to the reference sequence
                to verify variants and masked regions.''',
        choices=['mafft', 'edlib'])
    report_group.add_argument('--export_vcf',
        help='''Export a new vcf file named <sample>.final.vcf to the sample's results directory that
                contains all confirmed or introduced SNVs that were detected in the final consensus sequence.
                N-Masked variants are not written to the vcf file. No vcf file is exported for samples 
                missing a final consensus sequence.''',
        action='store_true')

    parser = add_help_group_to_parser(parser)
    args = parser.parse_args(remaining_argv)
    args = check_arguments(args)
    create_sample_reports(args)

if __name__ == '__main__':
    logger.error('This script is not intended to be directly executed with python. ' +\
                 'Please install package and use the package scripts instead.')
    exit(1)