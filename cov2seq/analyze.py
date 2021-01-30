import os, sys
import argparse
import json
import time
import pandas as pd
import numpy as np
import pdb
from clint.textui import colored, puts, indent

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

def scan_input_directories(input_directories, primer_schemes_dir, exclude, restrict,
                           min_len, max_len, min_qual, normalize, default_scheme, 
                           default_pore_model, scheme_specifics):
    samples = {}
    for input_dir in input_directories:
        if not os.path.exists(input_dir):
            logger.error("Input directory not found: {}".format(input_dir))
            exit(1)
        conf = os.path.join(input_dir, "run_configuration.json")
        if not os.path.exists(conf):
            logger.error("Run configuration file not found: {}".format(conf))
            exit(1)

        with open(conf, "r") as f:
            conf_data = json.load(f)

        run_id = conf_data['title']

        if 'poreModel' in conf_data:
            pore_model = conf_data['poreModel']
        else:
            logger.warning('Field "poreModel" missing in configuration file {}, assuming "{}"'.format(conf, default_pore_model))
            pore_model = default_pore_model

        primer_conf = os.path.join(input_dir, "primers.json")
        if not os.path.exists(primer_conf):
            scheme = default_scheme
            logger.warning('Primer configuration file not found: {}, assuming "{}"'.format(primer_conf, scheme))
        else:
            with open(primer_conf, "r") as f:
                primer_conf_data = json.load(f)
            scheme = primer_conf_data['name']
        scheme_dir = os.path.join(primer_schemes_dir, scheme)
        if not os.path.exists(scheme_dir):
            logger.error('Unknown primer scheme: "{}" not in {}'.format(scheme, primer_schemes_dir))
            exit(1)

        for sample_dict in conf_data['samples']:
            sample = sample_dict['name']
            barcodes = sample_dict['barcodes']
            if "ntc" in sample:
                continue
            if sample in exclude:
                continue
            if restrict and sample not in restrict:
                continue

            if not sample in samples:
                samples[sample] = {}
            if not scheme in samples[sample]:
                samples[sample][scheme] = {}
            samples[sample][scheme][run_id] = {
                'run_dir': input_dir,
                'barcodes': barcodes,
                'min_qual': scheme_specifics['min_qual'].get(scheme, min_qual),
                'min_len': scheme_specifics['min_len'].get(scheme, min_len),
                'max_len': scheme_specifics['max_len'].get(scheme, max_len),
                'pore_model': pore_model
                }
    return samples

def run_artic_medaka(args, scheme_specifics):
    retval = os.system("which artic medaka")
    if retval != 0:
        logger.error('Artic medaka command not found.')
        exit(1)

    cmds = []
    samples = scan_input_directories(args.input_directories, args.primer_schemes_dir, args.exclude, args.restrict,
                                     args.min_len, args.max_len, args.min_qual, args.normalize, args.scheme, 
                                     args.pore_model, scheme_specifics)
    if args.overview:
        data = []
        for sample in samples:
            for scheme in samples[sample]:
                for run_id in samples[sample][scheme]:
                    for barcode in samples[sample][scheme][run_id]['barcodes']:
                        data.append({
                            'sample': sample,
                            'scheme': scheme,
                            'run_id': run_id,
                            'barcode': barcode,
                            'min_qual': samples[sample][scheme][run_id]['min_qual'],
                            'min_len': samples[sample][scheme][run_id]['min_len'],
                            'max_len': samples[sample][scheme][run_id]['max_len'],
                            'pore_model': samples[sample][scheme][run_id]['pore_model'],
                            'run_dir': samples[sample][scheme][run_id]['run_dir']
                            })
        if args.index == 'sample':
            index = ['sample', 'run_id', 'barcode']
        else:
            index = ['run_id', 'sample', 'barcode']
        df = pd.DataFrame.from_records(data, index=index).sort_index()
        df.to_csv(args.overview, sep='\t')
        exit()

    for sample in samples:
        for scheme in samples[sample]:
            res_dir = os.path.abspath(os.path.join(args.nanopore_dir, sample, "artic-medaka", scheme)) + "/"
            if os.path.exists(res_dir):
                if not args.overwrite:
                    logger.error("Results directory {} already exists and --overwrite option is not set.".format(res_dir))
                    exit(1)
                cmd = "rm -r {}/*".format(res_dir)
                cmds.append(cmd)
            else:
                cmd = "mkdir -p {}".format(res_dir)
                cmds.append(cmd)

            joined_fq = os.path.abspath(os.path.join(res_dir, "{}.analyzed_reads.fastq").format(sample))
            pore_models = set()
            for run_id in samples[sample][scheme]:
                # copy primer configuration file to result dir
                primer_conf = os.path.join(samples[sample][scheme][run_id]['run_dir'], "primers.json")
                if os.path.exists(primer_conf):
                    cmd = "cp {} {}".format(primer_conf, res_dir)
                    cmds.append(cmd)

                # running artic guppyplex on each barcode dir and join resulting fastq files
                for bc in samples[sample][scheme][run_id]['barcodes']:
                    demux_dir = os.path.abspath(os.path.join(samples[sample][scheme][run_id]['run_dir'], "fastq_pass", bc)) + "/"
                    if not os.path.exists(demux_dir):
                        logger.error("Demux directory not found: {}".format(demux_dir))
                        exit(1)
                    run_fq = os.path.join(res_dir, "{}_q{}_{}-{}nt_{}.fastq".format(
                        run_id,
                        samples[sample][scheme][run_id]['min_qual'],
                        samples[sample][scheme][run_id]['min_len'],
                        samples[sample][scheme][run_id]['max_len'],
                        bc))
                    cmd = 'artic guppyplex --quality {} --min-length {} --max-length {} --directory {} --output {}'.format(
                        samples[sample][scheme][run_id]['min_qual'],
                        samples[sample][scheme][run_id]['min_len'],
                        samples[sample][scheme][run_id]['max_len'],
                        demux_dir,
                        run_fq)
                    cmds.append(cmd)

                    cmd = "cat {} >> {}".format(run_fq, joined_fq)
                    cmds.append(cmd)

                    pore_models.add(samples[sample][scheme][run_id]['pore_model'])

            # The new versions of medaka require a pore model. Therefore, sequencing experiments 
            # performed on different flowcell types cannot be analyzed together
            if len(pore_models) > 1:
                logger.error("Selected runs of sample {} where performed on flowcells with at least two different pore models: {}".format(sample, ", ".join(list(pore_models))))
                exit(1)

            # run artic minion medaka
            cmd = 'cd {} ; artic minion --medaka --medaka-model {} --normalise {} --threads {} --scheme-directory {} --read-file {} "{}" {}'.format(
                res_dir,
                pore_models.pop(),
                args.normalize,
                args.threads,
                args.primer_schemes_dir,
                joined_fq,
                scheme,
                sample)
            cmds.append(cmd)

            # writing artic version to file
            cmd = "artic --version > {}".format(os.path.join(res_dir, "artic.version.info"))
            cmds.append(cmd)

            # creating coverage files
            for suffix in ("sorted.bam", "trimmed.rg.sorted.bam", "primertrimmed.rg.sorted.bam"):
                cmd = 'bedtools genomecov -d -ibam {} > {}'.format(
                    os.path.join(res_dir, "{}.{}".format(sample, suffix)),
                    os.path.join(res_dir, "{}.{}.cov".format(sample, suffix)))
                cmds.append(cmd)

            # linking reference to result directory
            cmd = 'ln -s {}/{}/{}.reference.fasta {}/{}.reference.fasta'.format(
                args.primer_schemes_dir, scheme, scheme.split('/')[0], res_dir, scheme.split('/')[0])
            cmds.append(cmd)

            # chmod to add group write permissions
            dirname = os.path.abspath(os.path.join(args.nanopore_dir, sample))
            cmd = "chmod -R g+w {}".format(dirname)
            cmds.append(cmd)

    for cmd in cmds:
        print(colored.green("Running: ") + cmd, file=sys.stderr)
        if not args.dry_run:
            timerStart = time.perf_counter()
            retval = os.system(cmd)
            if retval != 0:
                print(colored.red('Command failed:' ) + cmd, file=sys.stderr)
                raise SystemExit(20)
            timerStop = time.perf_counter()