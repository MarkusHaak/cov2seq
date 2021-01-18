import os
import sys
from glob import glob
import json
from Bio import SeqIO
import numpy as np
import pandas as pd
import re
import subprocess
import vcf
from .verify import mafft, variants_from_alignment

import logging
logger = logging.getLogger()

def subdir_paths(basedir, level=2):
    basedir = basedir.rstrip(os.path.sep)
    return [d[len(basedir)+1:] for d in glob(os.path.join(*([basedir] + ['*'] * level))) if os.path.isdir(d)]

def load_reference(reference_fn, reference_annotation_fn):
    logger.info('Loading reference.')
    for record in SeqIO.parse(reference_fn, "fasta"):
        reference = record
        break
    with open(reference_annotation_fn, "r") as f:
        data = json.load(f)
    reference_genes = data['genes']
    return reference, reference_genes

def get_nanopore_runs(artic_dir):
    runs_info = []
    run_fastq_files = glob(os.path.join(artic_dir, '*_q*_*-*nt_*.fastq'))
    for fq in run_fastq_files:
        fn = os.path.basename(fq)
        m = re.search('(.+)_q(\d+)_(\d+)-(\d+)nt_(.*).fastq', fn)
        if not m:
            continue
        run, min_qual, min_len, max_len, barcode = m.group(1, 2, 3, 4, 5)
        qual, min_len, max_len = int(min_qual), int(min_len), int(max_len)
        runs_info.append((run, barcode, min_qual, min_len, max_len))
    return runs_info

def load_nanopore_info(samples, nanopore_dir):
    logger.info('Loading nanopore run information for selected samples.')
    nanopore_runs_data = []
    artic_runs_data = []
    for sample in samples:
        logger.debug('{}'.format(sample))
        basedir = os.path.join(nanopore_dir, sample, 'artic-medaka') #TODO: allow other than medaka
        sample_schemes = subdir_paths(basedir, level=2)
        if len(sample_schemes) == 0:
            logger.error("No nanopore results for sample {}".format(sample))
            exit(1)
        for scheme in sample_schemes:
            data_dir = os.path.join(basedir, scheme)
            # retrieve sequencing run information
            runs_info = get_nanopore_runs(data_dir)
            for run, barcode, min_qual, min_len, max_len in runs_info:
                nanopore_runs_data.append((sample, scheme, run, barcode, min_qual, min_len, max_len))

        # there must be a single artic directory that holds the merged results off all
        # nanopore runs
        if len(sample_schemes) > 1:
            data_dir = os.path.join(basedir, 'merged')
            if not os.path.isdir(data_dir):
                logger.error("Multiple different primer schemes for sample {}, but no merged results directory {} was found.".format(sample, data_dir))
                exit(1)
        else:
            data_dir = os.path.join(basedir, sample_schemes[0])
        # retrieve artic run information
        retval = subprocess.check_output(
            ["samtools","view","-c","-F","260", os.path.join(data_dir, sample+".sorted.bam")])
        mapped = int(retval.decode("utf-8").strip())
        retval = subprocess.check_output(
            ["samtools","view","-c","-F","260", os.path.join(data_dir, sample+".trimmed.rg.sorted.bam")])
        mapped_trimmed = int(retval.decode("utf-8").strip())
        retval = subprocess.check_output(
            ["samtools","view","-c","-F","260", os.path.join(data_dir, sample+".primertrimmed.rg.sorted.bam")])
        mapped_primertrimmed = int(retval.decode("utf-8").strip())
        artic_runs_data.append( (sample, data_dir, len(sample_schemes), len(runs_info), mapped, mapped_trimmed, mapped_primertrimmed) )

    nanopore_runs = pd.DataFrame(nanopore_runs_data, 
                                 columns=['sample', 'scheme', 'run', 'barcode', 'min_qual', 'min_len', 'max_len'])
    nanopore_runs = nanopore_runs.set_index('sample')
    artic_runs = pd.DataFrame(artic_runs_data,
                              columns=['sample', 'artic_dir', 'schemes', 'runs', 'mapped', 'mapped_trimmed', 'mapped_primertrimmed'])
    artic_runs = artic_runs.set_index('sample')
    assert len(artic_runs.index.duplicated() == 0)
    return artic_runs, nanopore_runs

def load_primer_schemes(primer_schemes_dir, primer_schemes):
    logger.info('Loading primer schemes')
    primers = {}
    amplicons = {}

    for scheme in primer_schemes:
        logger.debug('processing scheme {}'.format(scheme))
        scheme_name, scheme_version = scheme.split('/')

        scheme_tsv = os.path.join(primer_schemes_dir, scheme_name, scheme_version, scheme_name + ".tsv")
        scheme_bed = os.path.join(primer_schemes_dir, scheme_name, scheme_version, scheme_name + ".bed")
        if not os.path.exists(scheme_tsv):
            logger.error('.tsv file missing for primer scheme {}'.format(scheme))
            exit(1)
        if not os.path.exists(scheme_bed):
            logger.error('.tsv file missing for primer scheme {}'.format(scheme))
            exit(1)

        df1 = pd.read_csv(scheme_tsv, 
                          sep='\t', 
                          header=0, 
                          names=["name", "pool", "seq", "length", "%gc", "tm"], 
                          index_col="name", 
                          dtype={"pool":str, "seq":str, "length":int, "%gc":float, "tm":float})
        df2 = pd.read_csv(scheme_bed, 
                          sep='\t', 
                          header=None, 
                          names=["reference", "start", "end", "name", "pool", "strand"], 
                          index_col="name", 
                          dtype={"start":int, "end":int})
        primers[scheme] = pd.concat([df1, df2], axis=1)
        primers[scheme] = primers[scheme].loc[:,~primers[scheme].columns.duplicated()]
        primer_abrev = []
        amplicons_ = []
        for name,row in primers[scheme].iterrows():
            if "alt" in name.lower():
                s,num,ori,alt = name.split("_")
            else:
                s,num,ori = name.split("_")
                alt = ""
            primer_abrev.append(num+ori[0]+alt[3:])
            amplicons_.append(s + "_" + num)
        primers[scheme].insert(0,'abbrev', primer_abrev)
        primers[scheme]['amplicon'] = amplicons_
        
        amplicon_names = list(set(["_".join(n.split("_")[:2]) for n in list(primers[scheme].index)]))
        amplicon_names.sort(key=lambda x: int(x.split("_")[-1]))
        
        amplicon_data = []
        for amplicon in amplicon_names:
            indices = [name for name in list(primers[scheme].index) if name.startswith(amplicon+"_")]
            pstart = primers[scheme].loc[indices][primers[scheme].loc[indices]['strand'] == '+']['start'].min()
            start = primers[scheme].loc[indices][primers[scheme].loc[indices]['strand'] == '+']['end'].max() + 1
            end = primers[scheme].loc[indices][primers[scheme].loc[indices]['strand'] == '-']['start'].min() - 1
            pend = primers[scheme].loc[indices][primers[scheme].loc[indices]['strand'] == '-']['end'].max()
            pool = primers[scheme].loc[indices]['pool'].iloc[0]
            mingc = primers[scheme].loc[indices]['%gc'].min()
            amplicon_data.append( (amplicon, pool, pstart, start, end, pend, mingc) )
        amplicons[scheme] = pd.DataFrame(amplicon_data, 
                                         columns=['name', 'pool', 'pstart', 'start', 'end', 'pend', 'mingc'],
                                         index=amplicon_names)
    return primers, amplicons

def load_clades_info(nextstrain_ncov):
    logger.info('Loading nextstrain clades information.')
    clades_fn = os.path.join(nextstrain_ncov, "defaults", "clades.tsv")
    clades_df = pd.read_csv(clades_fn, sep="\t", header=0, skip_blank_lines=True).dropna()
    clades = list(clades_df['clade'].drop_duplicates())
    clades.sort()
    for i,clade in enumerate(clades):
        clades_df.loc[clades_df['clade'] == clade, 'color'] = "C{}".format(i+1)
    return clades_df

def cov_from_sorted_bam(bam_fn, cov_fn):
    if os.path.exists(cov_fn):
        if not os.access(bam_fn, os.W_OK):
            logger.error('Write permissions required for file {}'.format(cov_fn))
            exit(1)
    else:
        bam_dir = os.path.dirname(bam_fn)
        if not os.access(bam_dir, os.W_OK):
            logger.error('Write permissions required for directory {}'.format(bam_dir))
            exit(1)
    cmds = []
    cmd = 'bedtools genomecov -d -ibam {} > {}'.format(bam_fn, cov_fn)
    cmds.append(cmd)
    cmd = "chmod -R g+w {}".format(fn)
    cmds.append(cmd)
    for cmd in cmds:
        retval = os.system(cmd)
        if retval != 0:
            logger.error('Command failed:\n{}'.format(cmd))
            exit(1)

def get_coverage_from_bam(bam_fn):
    cov_fn = bam_fn + '.cov'
    if not os.path.exists(bam_fn):
        cov_fn = cov_from_sorted_bam(bam_fn, cov_fn)
    return np.genfromtxt(cov_fn, delimiter='\t', usecols=(2), dtype=np.int32)

def get_nanopore_coverage(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    bam_fn = os.path.join(artic_dir, "{}.sorted.bam.cov".format(sample))
    return get_coverage_from_bam(bam_fn)

def get_nanopore_coverage_trimmed(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']   
    bam_fn = os.path.join(artic_dir, "{}.trimmed.rg.sorted.bam.cov".format(sample))
    return get_coverage_from_bam(bam_fn)

def get_nanopore_coverage_primertrimmed(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    bam_fn = os.path.join(artic_dir, "{}.primertrimmed.rg.sorted.bam.cov".format(sample))
    return get_coverage_from_bam(bam_fn)

def get_nanopore_pool_coverage(sample, artic_runs, nanopore_runs, amplicons, reference):
    sample_schemes = nanopore_runs.loc[sample, 'scheme']
    if type(sample_schemes) == str:
        sample_schemes = [sample_schemes]
    else:
        sample_schemes = list(sample_schemes.drop_duplicates())
    coverage_pool = {scheme: {} for scheme in sample_schemes}
    for scheme in sample_schemes:
        for pool in list(amplicons[scheme]['pool'].drop_duplicates()):

            bam_fn = os.path.join(artic_dir, "{}.primertrimmed.{}.sorted.bam".format(sample, pool))
            if os.path.exists(bam_fn):
                coverage_pool[scheme][pool] = get_coverage_from_bam(bam_fn)
            else:
                logger.warning(('Bam file {} missing. Coverage of primer scheme {} pool {} ' + \
                                'will be displayed as zero which might not be the case.').format(bam_fn, scheme, pool))
                coverage_pool[scheme][pool] = np.zeros(shape=(len(reference.seq),), dtype=np.int32)
    return coverage_pool

def get_nanopore_coverage_mask(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    fn = os.path.join(os.path.join(artic_dir, "{}.coverage_mask.txt".format(sample)))
    cov_mask = np.genfromtxt(fn, delimiter='\t', usecols=(1,2), dtype=np.int32)
    if len(cov_mask.shape) == 1:
        cov_mask = cov_mask[None, :]
    return cov_mask

def get_illumina_coverage_and_mappings(sample, illumina_dir, reference):
    coverage_illumina = np.zeros(shape=(len(reference.seq),), dtype=np.int32)
    mapped_illumina = []
    sample_ill_dir = os.path.join(illumina_dir, sample)
    if os.path.exists(sample_ill_dir):
        for bam_fn in glob(os.path.join(sample_ill_dir, '*.bam')):
            if 'extended' in bam_fn or 'sorted' in bam_fn:
                continue
            cov_fn = bam_fn + '.cov'
            bam_sorted_fn = bam_fn[:-3] + 'sorted.bam'
            if not os.path.exists(cov_fn):
                if not os.access(sample_ill_dir, os.W_OK):
                    logger.warning('''Unable to account illumina coverage of {},
                                      no write permissions for {}'''.format(bam_fn, sample_ill_dir))
                    continue
                if not os.path.exists(bam_sorted_fn):
                    cmd = "samtools sort {} > {}".format(bam_fn, bam_sorted_fn)
                    os.system(cmd)
                    cmd = "chmod g+x {}".format(bam_sorted_fn)
                    os.system(cmd)
                cmd = "bedtools genomecov -d -ibam {} > {}".format(bam_sorted_fn, cov_fn)
                os.system(cmd)
                cmd = "chmod g+x {}".format(cov_fn)
                os.system(cmd)
            coverage = np.genfromtxt(cov_fn, delimiter='\t', usecols=(2), dtype=np.int32)
            if coverage.shape != coverage_illumina.shape:
                logger.warning('''Shapes of coverage extracted from {} and reference do not match
                                  possibly indicating that a different reference was used 
                                  during mapping of Illumina reads.'''.format(cov_fn))
            coverage_illumina += coverage
            # mappings
            retval = subprocess.check_output("samtools view -c -F 260 {}".format(bam_sorted_fn).split(" "))
            mapped_illumina.append((os.path.basename(cov_fn), 
                                    int(retval.decode("utf-8").strip()), 
                                    np.sum(coverage)))
    mapped_illumina = pd.DataFrame(mapped_illumina, 
                                   columns=['bam_file', 'mapped_reads', 'mapped_bases'])
    mapped_illumina = mapped_illumina.set_index('bam_file', keep=False)
    return coverage_illumina, mapped_illumina

def approximate_sanger_coverage(sample, sanger_dir, reference, amplicons, primers):
    coverage_sanger = np.zeros(shape=(len(reference.seq),), dtype=np.int32)
    sample_sanger_dir = os.path.join(sanger_dir, sample)
    if os.path.exists(sample_sanger_dir):
        for ab1_fn in glob(os.path.join(sample_sanger_dir, '*.ab1')):
            m = re.search('.*-(.*).ab1', os.path.basename(ab1_fn))
            if not m:
                logger.warning('''Sanger trace file has wrong naming scheme
                                  and could therefore not be parsed: {}'''.format(ab1_fn))
                continue
            abbrev = m.group(1)
            amplicon = primers['nCoV-2019/V3'].loc[(primers['nCoV-2019/V3']['abbrev'] == abbrev), 'amplicon'].to_list()[0]
            start, end = amplicons['nCoV-2019/V3'].loc[amplicon, ['start', 'end']].to_list()
            coverage_sanger[start:end] += 1
    return coverage_sanger

def assign_clade(sample, artic_runs, results_dir, nextstrain_ncov, repeat_assignment=True):
    '''performs a clade assignment using nextstrain's assign_clade.py script and returns
    the results as a tuple of three values: 'fasta header', 'clade' and 'parent clade' '''
    consensus_fn = os.path.join(results_dir, sample, '{}.final.fasta'.format(sample))
    clade_fn = os.path.join(results_dir, sample, '{}.clade.tsv'.format(sample))
    version_fn = os.path.join(results_dir, sample, 'assign_clades.py.version.info')
    if not os.path.exists(consensus_fn):
        artic_dir = artic_runs.loc[sample, 'artic_dir']
        consensus_fn = os.path.join(artic_dir, '{}.consensus.fasta'.format(sample))
        clade_fn = os.path.join(artic_dir, '{}.clade.tsv'.format(sample))
        version_fn = os.path.join(artic_dir, 'assign_clades.py.version.info')
    # repeat clade assignment by default ?
    if not os.path.exists(clade_fn) or repeat_assignment:
        cmd =  'cd {} && '.format(nextstrain_ncov) 
        cmd += '{} scripts/assign_clades.py --sequences {} --output {} && '.format(sys.executable, consensus_fn, clade_fn)
        cmd += 'git log --pretty="%H" -n 1  > {}'.format(version_fn)
        ret = os.system(cmd)
        if ret != 0:
            logger.error("Clade assignment command failed: {}".format(cmd))
            exit(1)
        cmd = 'chmod g+w {}'.format(clade_fn)
        os.system(cmd)
        cmd = 'chmod g+w {}'.format(version_fn)
        os.system(cmd)
    try:
        df = pd.read_csv(clade_fn, header=0, sep="\t")
        if len(df) == 0:
            return ('', 'unknown', '')
        return tuple(df.iloc[0])
    except:
        logger.warning("Clade assignment file is of unknown format: {}".format(clade_fn))
        return ('', 'unknown', '')

def run_extended_snv_pipeline(sample, artic_runs, reference_fn, reference_annotation_fn):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    sorted_bam_fn = os.path.join(artic_dir, '{}.primertrimmed.rg.sorted.bam'.format(sample))
    merged_vcf_fn = os.path.join(artic_dir, '{}.merged.vcf.gz'.format(sample))
    annotated_vcf_fn = os.path.join(artic_dir, '{}.merged.ann.vcf'.format(sample))
    strand_bias_vcf_fn = os.path.join(artic_dir, "{}.longshot.01.vcf".format(sample))
    for fn in [annotated_vcf_fn, strand_bias_vcf_fn]:
        if os.path.exists(fn):
            if not os.access(fn, os.W_OK):
                logger.warning('Write permissions required for file {}'.format(fn))
                return False
    cmds = []
    # re-run longshot variant calling with strand bias filter
    cmd = 'longshot -P 0.01 -F -A --no_haps --bam {} --ref {} --out {} --potential_variants {}'.format(
        sorted_bam_fn, reference_fn, strand_bias_vcf_fn, merged_vcf_fn)
    cmds.append(cmd)
    # run annotation of merged vcf
    cmd = 'zcat {0} > {1} ; vcf-annotator {1} {2} --output {1}'.format(
        merged_vcf_fn, annotated_vcf_fn, reference_annotation_fn)
    cmds.append(cmd)
    # chmod
    for fn in [annotated_vcf_fn, strand_bias_vcf_fn]:
        if oct(os.stat(fn).st_mode)[-2] not in '732':
            cmd = "chmod -R g+w {}".format(fn)
            cmds.append(cmd)
    for cmd in cmds:
        retval = os.system(cmd)
        if retval != 0:
            logger.error('Command failed:\n{}'.format(cmd))
            exit(1)
    return True

def parse_vcf(vcf_fn, info_fcts={}, constant_fields={}):
    vcf_reader = vcf.Reader(filename=vcf_fn)
    vcf_data = []
    index = []
    for record in vcf_reader:
        fields = {}
        # multiple genotypes not supported
        if len(record.ALT) > 1:
            logger.info("Skipping VCF {} {} (annotated vcf: multiple genotypes not supported)".format(sample, record))
            continue
        fields['site'], fields['ref'], fields['alt'] = int(record.POS), str(record.REF), str(record.ALT[0])
        fields['qual'], fields['filter'] = record.QUAL, ",".join(record.FILTER)
        for key, f in info_fcts.items():
            fields[key] = f(record.INFO)
        fields.update(constant_fields)
        
        index.append("{}{}{}".format(fields['ref'], fields['site'], fields['alt']))
        vcf_data.append(fields)
    return pd.DataFrame(vcf_data, index=pd.Index(index))

def load_snv_info(sample, artic_runs, results_dir, reference_fn, reference_annotation_fn, clades_df):
    snv_info = []
    multiindex_rows = [[],[]]
    annotator_fcts = {
        'Pool' : lambda info: info['Pool'],
        'AminoAcidChange' : lambda info: ", ".join(str(e) if e else '' for e in info['AminoAcidChange']),
        'RefCodon' : lambda info: " ".join(str(e) if e else '' for e in info['RefCodon']),
        'AltCodon' : lambda info: " ".join(str(e) if e else '' for e in info['AltCodon']),
        'Gene' : lambda info: ", ".join(str(e) if e else '' for e in info['Gene']),
        'Product' : lambda info: ", ".join(str(e) if e else '' for e in info['Product'])
    }
    longshot_fcts = {
        'cov' : lambda info: info['AC'][0] + info['AC'][1] + info['AM'],
        '#ref' : lambda info: info['AC'][0],
        '#alt' : lambda info: info['AC'][1],
        '#amb' : lambda info: info['AM']
    }
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    annotated_vcf_fn = os.path.join(artic_dir, "{}.merged.ann.vcf".format(sample))
    pass_vcf_fn = os.path.join(artic_dir, "{}.pass.vcf.gz".format(sample))
    fail_vcf_fn = os.path.join(artic_dir, "{}.fail.vcf".format(sample))
    longshot_vcf_fn = os.path.join(artic_dir, "{}.longshot.vcf".format(sample))
    strand_bias_vcf_fn = os.path.join(artic_dir, "{}.longshot.01.vcf".format(sample))
    sample_final_consensus = os.path.join(results_dir, sample, "{}.final.fasta".format(sample))
    # perform snv analysis if vcf files are missing
    if (not os.path.exists(annotated_vcf_fn)) or (not os.path.exists(strand_bias_vcf_fn)):
        if not run_extended_snv_pipeline(sample, artic_runs, reference_fn, reference_annotation_fn):
            logger.warning('Failed to run extended SNV analysis for sample {}'.format(sample))
    # parse vcf files
    vcf_ann = parse_vcf(annotated_vcf_fn, info_fcts=annotator_fcts)
    vcf_pass = parse_vcf(pass_vcf_fn, constant_fields={'snv_filter': True}, info_fcts=longshot_fcts)
    vcf_fail = parse_vcf(fail_vcf_fn, constant_fields={'snv_filter': False}, info_fcts=longshot_fcts)
    vcf_longshot = parse_vcf(longshot_vcf_fn, info_fcts=longshot_fcts)
    vcf_bias = parse_vcf(strand_bias_vcf_fn, info_fcts=longshot_fcts)
    if len(vcf_pass) and len(vcf_fail):
        for index in vcf_pass.index.intersection(vcf_fail):
            logger.warning('SNV {} is both in the pass and fail .vcf files created by the ARTIC snv_filter tool')
    # extract medaka variant information
    vcf_medaka = vcf_ann[['Pool', 'site', 'ref', 'alt', 'qual', 'filter']]
    vcf_medaka.columns = pd.MultiIndex.from_product([['medaka variant'], vcf_medaka.columns])
    # extract and format vcf-annotator information
    vcf_ann["Product"] = vcf_ann["Product"].str.replace('[space]', ' ', regex=False)
    vcf_annotation = vcf_ann[['AminoAcidChange', 'RefCodon', 'AltCodon', 'Gene', 'Product']]
    vcf_annotation.columns = pd.MultiIndex.from_product([['vcf-annotator'], vcf_annotation.columns])
    vcf_annotation = vcf_annotation.drop_duplicates()
    # extract information about which SNVs passed and failed the ARTIC vcf_filter
    vcf_artic_filter = pd.concat([vcf_pass,vcf_fail])['snv_filter'].to_frame()
    vcf_artic_filter.columns = pd.MultiIndex.from_product([['ARTIC'], vcf_artic_filter.columns])
    vcf_artic_filter = vcf_artic_filter.reset_index()
    vcf_artic_filter = vcf_artic_filter.drop_duplicates().set_index('index')
    if np.any(vcf_artic_filter.index.duplicated()):
        duplicated = list(vcf_artic_filter.index[vcf_artic_filter.index.duplicated()])
        logger.error('Variant(s) {} of sample {} are classified as both pass and fail by ARTIC snv_filter'.format(duplicated, sample))
        exit(1)
    # extract information from longshot and decorate it with information about the longshot strand bias filter
    common_columns = list(vcf_longshot.columns.values)
    common_columns.append('index')
    vcf_bias['strand bias'] = False
    vcf_longshot_bias = pd.merge(vcf_longshot.reset_index(), 
                                 vcf_bias.reset_index(), 
                                 how='left', 
                                 left_on=common_columns, 
                                 right_on=common_columns).fillna(True).set_index('index')
    vcf_longshot_bias = vcf_longshot_bias[['qual', 'filter', 'cov', '#ref', '#alt', '#amb', 'strand bias']]
    vcf_longshot_bias.columns = pd.MultiIndex.from_product([['longshot'], vcf_longshot_bias.columns])
    # if a final consensus sequence is present inthe results directory,
    # align it against the reference to extract information about SNVs that were confirmed or rejected
    if os.path.exists(sample_final_consensus):
        alignment = mafft([sample_final_consensus], reference_fn)
        vcf_confirmed = variants_from_alignment(alignment, [sample]).loc[sample]
        vcf_confirmed['confirmation'] = True
        vcf_confirmed = vcf_confirmed['confirmation'].to_frame()
        vcf_confirmed.columns = pd.MultiIndex.from_product([['final'], vcf_confirmed.columns])
    else:
        vcf_confirmed = pd.DataFrame([], columns=pd.MultiIndex.from_arrays([['final'],['confirmation']]))

    # merge information from individual vcf files to one table
    if not vcf_medaka.loc[vcf_longshot_bias.index.drop_duplicates()].index.equals(vcf_longshot_bias.index):
        logger.warning('Number of variants of potentially different pools identified by ' + \
                       'medaka variant does not match the number of entries ' + \
                       'in the longshot output. Medaka variant information is therefore ' + \
                       'dropped for variants where an unambiguous assignment is not possible.')
        reconstruction = []
        for index in vcf_longshot_bias.index.drop_duplicates():
            df_medaka = vcf_medaka.loc[index]
            if type(df_medaka) == pd.core.series.Series:
                df_medaka = df_medaka.to_frame().T
            df_longshot_bias = vcf_longshot_bias.loc[index]
            if type(df_longshot_bias) == pd.core.series.Series:
                df_longshot_bias = df_longshot_bias.to_frame().T
            if len(df_medaka) == len(df_longshot_bias):
                reconstruction.append(pd.concat([df_medaka, df_longshot_bias], axis=1))
            else:
                df_medaka = df_medaka.iloc[:len(df_longshot_bias)]
                df_medaka[('medaka variant', 'Pool')] = 'N/A'
                df_medaka[('medaka variant', 'qual')] = np.nan
                reconstruction.append(pd.concat([df_medaka, df_longshot_bias], axis=1))
        snv_info = pd.concat(reconstruction, axis=0, sort=False)
    else:
        snv_info = pd.concat([vcf_medaka.loc[vcf_longshot_bias.index.drop_duplicates()], vcf_longshot_bias], axis=1)
    columns = snv_info.columns #save column order, since pd.concat sorts the column names
    snv_info = pd.concat([snv_info, vcf_medaka.loc[vcf_medaka.index.difference(vcf_longshot_bias.index.drop_duplicates())]]).sort_values(('medaka variant','site'))
    snv_info = snv_info.loc[:, columns] #restore column order
    snv_info = snv_info.join(vcf_artic_filter, how='left')
    snv_info = snv_info.join(vcf_annotation, how='left')
    #if os.path.exists(sample_final_consensus):
    #    breakpoint()
    snv_info = snv_info.join(vcf_confirmed, how='outer')
    sel_null = pd.isnull(snv_info[('medaka variant', 'site')])
    if np.any(sel_null):
        snv_info.loc[sel_null, ('medaka variant', 'site')] = int(snv_info.loc[sel_null].index.str.extract(r"\D+(\d+)\D+")[0])
        snv_info.loc[sel_null, ('medaka variant', 'ref')] = snv_info.loc[sel_null].index.str.extract(r"(\D+)\d+\D+")[0]
        snv_info.loc[sel_null, ('medaka variant', 'alt')] = snv_info.loc[sel_null].index.str.extract(r"\D+\d+(\D+)")[0]
    return snv_info.sort_values(('medaka variant', 'site'))

def get_software_versions(sample, artic_runs, results_dir):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    version_files = [
        ('artic fieldbioinformatics', os.path.join(artic_dir, 'artic.version.info'))]
    if os.path.exists(os.path.join(results_dir, sample, 'assign_clades.py.version.info')):
        version_files.append(('nextstrain assign_clades.py', 
                              os.path.join(results_dir, sample, 'assign_clades.py.version.info')))
    else:
        version_files.append(('nextstrain assign_clades.py', 
                              os.path.join(artic_dir, 'assign_clades.py.version.info')))
    software_versions = {}
    for software, fn in version_files:
        if os.path.exists(fn):
            with open(fn, 'r') as f:
                software_versions[software] = f.readline().strip()
        else:
            software_versions[software] = 'unknown'
    return software_versions
