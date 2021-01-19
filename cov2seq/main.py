import os
import sys
import configparser
import pathlib
import pkg_resources
import argparse
from .report import create_sample_reports

import logging
logger = logging.getLogger()
shandler = logging.StreamHandler()
formatter = logging.Formatter('%(levelname)s: %(message)s')
shandler.setFormatter(formatter)
logger.addHandler(shandler)
logger.setLevel(logging.DEBUG)

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

def read_configuration(argv, pkg_dir, defaults=None, fields=None):
    def update_defaults(conf_file):
        config = configparser.SafeConfigParser()
        config.read([conf_file])
        data_fields = ['DEFAULT']
        if fields:
            data_fields += fields
        for field in data_fields:
            defaults.update(dict(config.items(field)))

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
        update_defaults(args.config_file)

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
    main_group.add_argument('--reference_fn',
        help='Fasta file containing the reference sequence used in the artic pipeline.')
    main_group.add_argument('--reference_annotation_fn',
        help='.json file with basic gene annotation information.')
    main_group.add_argument('--primer_schemes_dir',
        help='''Direcotry containing the primer scheme files in a directory tree as follows:
             <primer_schemes_dir>/<scheme_name>/<scheme_version>/<scheme_name>.*''')
    main_group.add_argument('--results_dir',
        help='''Directory for saving sample results.''')
    main_group.add_argument('--nextstrain_ncov',
        help='''Path to the directory containing the Nextstrain ncov git repository. This should
             be updated regularly for latest clade assignments.''')
    return parser

def check_arguments(args):
    args.nanopore_dir = os.path.expanduser(args.nanopore_dir)
    args.primer_schemes_dir = os.path.expanduser(args.primer_schemes_dir)
    args.results_dir = os.path.expanduser(args.results_dir)
    args.nextstrain_ncov = os.path.expanduser(args.nextstrain_ncov)
    args.reference_fn = os.path.expanduser(args.reference_fn)
    args.reference_annotation_fn = os.path.expanduser(args.reference_annotation_fn)
    args.illumina_dir = os.path.expanduser(args.illumina_dir)
    args.sanger_dir = os.path.expanduser(args.sanger_dir)
    for arg in [args.nanopore_dir,
                args.primer_schemes_dir,
                args.results_dir,
                args.nextstrain_ncov]:
        if not os.path.isdir(arg):
            logger.error('Directory {} does not exist but is required.'.format(arg))
            exit(1)
        if not os.access(arg, os.W_OK):
            logger.error('Write permissions required for directory {}'.format(arg))
            exit(1)
    for arg in [args.reference_fn,
                args.reference_annotation_fn]:
        if not os.path.isfile(arg):
            logger.error('File {} does not exist but is required.'.format(arg))
            exit(1)
    for arg in [args.illumina_dir,
                args.sanger_dir]:
        if not os.path.isdir(arg):
            logger.warning('Directory {} does not exist.'.format(arg))
        else:
            if not os.access(arg, os.W_OK):
                logger.error('Write permissions required for directory {}'.format(arg))
                exit(1)
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

def init_parser(argv, defaults, script_descr="", fields=[]):
    pkg_dir, pkg_version, pkg_descr = get_package_info()
    
    nextstrain_ncov = os.path.join(pkg_dir, 'ncov')
    if os.path.exists(nextstrain_ncov):
        defaults.update({'nextstrain_ncov' : nextstrain_ncov})

    conf_parser, remaining_argv, defaults = read_configuration(argv, pkg_dir, defaults, fields)

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

def update(argv=None):
    if argv is None:
        argv = sys.argv
    defaults = {} # configuration file independent default values; lowest priority
    parser, remaining_argv = init_parser(argv, defaults)


    parser = add_help_group_to_parser(parser)
    args = parser.parse_args(remaining_argv)

def report(argv=None):
    if argv is None:
        argv = sys.argv
    defaults = {} # configuration file independent default values; lowest priority
    script_descr = "This script creates comprehensive reports for one or several samples."
    parser, remaining_argv = init_parser(argv, defaults, script_descr=script_descr, fields=['REPORT'])

    report_group = parser.add_argument_group('Report Option Group')
    report_group.add_argument('-s', '--samples',
        help='Relative glob patterns matching sample names present in <nanopore_dir> directory.',
        nargs='+',
        required=True)
    report_group.add_argument('--repeat_assignment',
        help='Repeat the clade assignment even if a clade assignment was already once performed.',
        action='store_true')
    report_group.add_argument('--threshold_limit',
        help='''Minimal acceptable nanopore coverage. Regions with a coverage below this 
                this threshold are highlighted red.''',
        type=int)
    report_group.add_argument('--threshold_low',
        help='''Minimal acceptable nanopore coverage to assume that variant calling based on nanopore 
                data alone works sufficiently well. Regions with a coverage below this 
                this threshold are highlighted orange.''',
        type=int)

    parser = add_help_group_to_parser(parser)
    args = parser.parse_args(remaining_argv)
    args = check_arguments(args)

    pkg_dir,_,_ = get_package_info()
    create_sample_reports(args, pkg_dir)

if __name__ == '__main__':
    print('''This script is not intended for standalone command-line use. 
             Please install package and use cov2seq-update or cov2seq-report instead.''')
    exit(1)