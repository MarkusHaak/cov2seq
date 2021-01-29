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
from tqdm import tqdm
from .verify import compare_consensus_to_reference

import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

aa_single_letter_ = {"Gly":"G", "Pro":"P",
                     "Ala":"A", "Val":"V",
                     "Leu":"L", "Ile":"I",
                     "Met":"M", "Cys":"C",
                     "Phe":"F", "Tyr":"Y",
                     "Trp":"W", "His":"H",
                     "Lys":"K", "Arg":"R",
                     "Gln":"Q", "Asn":"N",
                     "Glu":"E", "Asp":"D",
                     "Ser":"S", "Thr":"T"}

def selected_samples_from_patterns(nanopore_dir, patterns):
    sample_paths = []
    for pattern in patterns:
        sample_paths.extend(glob(os.path.join(nanopore_dir, pattern)))
    selected_samples = [os.path.basename(p)  for p in set(sample_paths)]
    selected_samples.sort()
    return selected_samples

def subdir_paths(basedir, level=2):
    basedir = basedir.rstrip(os.path.sep)
    return [d[len(basedir)+1:] for d in glob(os.path.join(*([basedir] + ['*'] * level))) if os.path.isdir(d)]

def load_reference(reference_fn):
    logger.info('Loading reference.')
    reference = next(SeqIO.parse(reference_fn, "gb"))
    data = []
    for feature in reference.features:
        if feature.type=="CDS":
            data.append({"gene" : feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else "",
                         "strand" : feature.strand,
                         "start" : feature.location.start.position,
                         "end" : feature.location.end.position,
                         "protein_id" : feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else "",
                         "product" : feature.qualifiers['product'][0] if 'product' in feature.qualifiers else "",
                         "note" : feature.qualifiers['note'][0] if 'note' in feature.qualifiers else ""})
    reference_genes = pd.DataFrame(data)
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

def load_nanopore_info(samples, nanopore_dir, sorted_reads=True, trimmed_reads=True, primertrimmed_reads=True):
    logger.info('Loading nanopore run information for selected samples.')
    nanopore_runs_data = []
    artic_runs_data = []
    for sample in tqdm(samples):
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
                logger.error("Multiple different primer schemes for sample {}, " +\
                             "but no merged results directory {} was found.".format(sample, data_dir))
                exit(1)
        else:
            data_dir = os.path.join(basedir, sample_schemes[0])
        sample_artic_dict = {'sample':sample,
                             'artic_dir':data_dir,
                             'schemes':len(sample_schemes),
                             'runs':len(nanopore_runs_data)}
        # retrieve artic run information
        if sorted_reads:
            retval = subprocess.check_output(
                ["samtools","view","-c","-F","260", os.path.join(data_dir, sample+".sorted.bam")])
            sample_artic_dict['mapped'] = int(retval.decode("utf-8").strip())
        if trimmed_reads:
            retval = subprocess.check_output(
                ["samtools","view","-c","-F","260", os.path.join(data_dir, sample+".trimmed.rg.sorted.bam")])
            sample_artic_dict['mapped_trimmed'] = int(retval.decode("utf-8").strip())
        if primertrimmed_reads:
            retval = subprocess.check_output(
                ["samtools","view","-c","-F","260", os.path.join(data_dir, sample+".primertrimmed.rg.sorted.bam")])
            sample_artic_dict['mapped_primertrimmed'] = int(retval.decode("utf-8").strip())
        artic_runs_data.append(sample_artic_dict)

    nanopore_runs = pd.DataFrame(nanopore_runs_data, 
                                 columns=['sample', 'scheme', 'run', 'barcode', 'min_qual', 'min_len', 'max_len'])
    nanopore_runs = nanopore_runs.set_index('sample')
    artic_runs = pd.DataFrame(artic_runs_data)
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

        scheme_insert_fn = os.path.join(primer_schemes_dir, scheme_name, scheme_version, scheme_name + ".insert.bed")
        scheme_primer_fn = os.path.join(primer_schemes_dir, scheme_name, scheme_version, scheme_name + ".primer.bed")

        if not os.path.exists(scheme_insert_fn):
            logger.error('File {} missing for primer scheme {}'.format(scheme_insert_fn, scheme))
            exit(1)
        if not os.path.exists(scheme_primer_fn):
            logger.error('File {} missing for primer scheme {}'.format(scheme_primer_fn, scheme))
            exit(1)
        
        # load primer data
        df_p = pd.read_csv(scheme_primer_fn, sep='\t', header=None,
                           names=["reference", "pstart", "pend", "name", "pool", "strand"])
        df_p['amplicon_number'] = df_p['name'].str.extract(r"{}_(\d+)_*".format(scheme_name)).astype(np.int32)
        df_p['abbrev'] = df_p['name'].str.extract(r"{}_(\d+_.*)".format(scheme_name))
        df_p['amplicon'] = scheme_name + '_' + df_p['amplicon_number'].astype(str)
        # load amplicon data
        df_i = pd.read_csv(scheme_insert_fn, sep='\t', header=None,
                           names=["reference", "start", "end", "number", "pool", "strand"])
        df_i['name'] = scheme_name + '_' + df_i['number'].astype(str)
        df_i = df_i.set_index('name')
        pstart = df_p.loc[df_p['strand'] == '+'].groupby('amplicon').apply(lambda x: x.loc[x.pstart.idxmin(), 'pstart'])
        pend = df_p.loc[df_p['strand'] == '-'].groupby('amplicon').apply(lambda x: x.loc[x.pend.idxmax(), 'pend'])
        pstart_primer = df_p.loc[df_p['strand'] == '+'].groupby('amplicon').apply(lambda x: x.loc[x.pstart.idxmin(), 'name'])
        pend_primer = df_p.loc[df_p['strand'] == '-'].groupby('amplicon').apply(lambda x: x.loc[x.pend.idxmax(), 'name'])
        df_i.loc[pstart.index, 'pstart'] = pstart
        df_i.loc[pstart.index, 'pend'] = pend
        df_i.loc[pstart.index, 'leftmost_primer'] = pstart_primer
        df_i.loc[pstart.index, 'rightmost_primer'] = pend_primer
        
        df_p = df_p.set_index('name', drop=True)
        primers[scheme] = df_p
        amplicons[scheme] = df_i
    return primers, amplicons

def load_clades_info(nextstrain_ncov_dir):
    logger.info('Loading nextstrain clades information.')
    clades_fn = os.path.join(nextstrain_ncov_dir, "defaults", "clades.tsv")
    subclades_fn = os.path.join(nextstrain_ncov_dir, "defaults", "subclades.tsv")
    clades_df = pd.read_csv(clades_fn, sep="\t", header=0, skip_blank_lines=True).dropna()
    subclades_df = pd.read_csv(clades_fn, sep="\t", header=None, names=list(clades_df.columns), skip_blank_lines=True).dropna()
    #clades = list(clades_df['clade'].drop_duplicates())
    #clades.sort()
    #for i,clade in enumerate(clades):
    #    clades_df.loc[clades_df['clade'] == clade, 'color'] = "C{}".format(i+1)
    return clades_df, subclades_df

def cov_from_sorted_bam(bam_fn, cov_fn):
    logger.debug("Extracting coverage from {} using bedtools".format(bam_fn))
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
    cmd = "chmod -R g+w {}".format(cov_fn)
    cmds.append(cmd)
    for cmd in cmds:
        retval = os.system(cmd)
        if retval != 0:
            logger.error('Command failed:\n{}'.format(cmd))
            exit(1)

def get_coverage_from_bam(bam_fn):
    cov_fn = bam_fn + '.cov'
    if not os.path.exists(cov_fn):
        cov_from_sorted_bam(bam_fn, cov_fn)
    return np.genfromtxt(cov_fn, delimiter='\t', usecols=(2), dtype=np.int32)

def get_nanopore_coverage(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    bam_fn = os.path.join(artic_dir, "{}.sorted.bam".format(sample))
    return get_coverage_from_bam(bam_fn)

def get_nanopore_coverage_trimmed(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']   
    bam_fn = os.path.join(artic_dir, "{}.trimmed.rg.sorted.bam".format(sample))
    return get_coverage_from_bam(bam_fn)

def get_nanopore_coverage_primertrimmed(sample, artic_runs):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    bam_fn = os.path.join(artic_dir, "{}.primertrimmed.rg.sorted.bam".format(sample))
    return get_coverage_from_bam(bam_fn)

def get_nanopore_pool_coverage(sample, artic_runs, nanopore_runs, amplicons, reference):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    sample_schemes = nanopore_runs.loc[sample, 'scheme']
    if type(sample_schemes) == str:
        sample_schemes = [sample_schemes]
    else:
        sample_schemes = list(sample_schemes.drop_duplicates())
    coverage_pools = {scheme: {} for scheme in sample_schemes}
    for scheme in sample_schemes:
        for pool in list(amplicons[scheme]['pool'].drop_duplicates()):

            bam_fn = os.path.join(artic_dir, "{}.primertrimmed.{}.sorted.bam".format(sample, pool))
            if os.path.exists(bam_fn):
                coverage_pools[scheme][pool] = get_coverage_from_bam(bam_fn)
            else:
                logger.warning(('Bam file {} missing. Coverage of primer scheme {} pool {} ' + \
                                'will be displayed as zero which might not be the case.').format(bam_fn, scheme, pool))
                coverage_pools[scheme][pool] = np.zeros(shape=(len(reference.seq),), dtype=np.int32)
    return coverage_pools

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
    mapped_illumina = mapped_illumina.set_index('bam_file', append=False)
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

def assign_clade(sample, artic_runs, results_dir, nextstrain_ncov_dir, repeat_assignment=True):
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
        cmd =  'cd {} && '.format(nextstrain_ncov_dir) 
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
        cmd = 'chmod g+w {}'.format(os.path.join(nextstrain_ncov_dir, 'clade_assignment_tmp*'))
        os.system(cmd)
    try:
        df = pd.read_csv(clade_fn, header=0, sep="\t")
        if len(df) == 0:
            return ('', 'unknown', '')
        return tuple(df.iloc[0])
    except:
        logger.warning("Clade assignment file is of unknown format: {}".format(clade_fn))
        return ('', 'unknown', '')

def run_extended_snv_pipeline(sample, artic_runs, snpeff_dir, reference_fasta_fn):
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    sorted_bam_fn = os.path.join(artic_dir, '{}.primertrimmed.rg.sorted.bam'.format(sample))
    merged_vcf_fn = os.path.join(artic_dir, '{}.merged.vcf.gz'.format(sample))
    annotated_vcf_fn = os.path.join(artic_dir, '{}.merged.snpeff.vcf'.format(sample))
    strand_bias_vcf_fn = os.path.join(artic_dir, "{}.longshot.01.vcf".format(sample))
    if os.path.exists(annotated_vcf_fn):
        if not os.access(annotated_vcf_fn, os.W_OK):
            logger.warning('Write permissions required for file {}'.format(annotated_vcf_fn))
            return False
    cmds = []
    if not os.path.exists(strand_bias_vcf_fn):
        # re-run longshot variant calling with strand bias filter
        if not os.path.exists(reference_fasta_fn + '.fai'):
            cmd = "samtools faidx {}".format(reference_fasta_fn)
            cmds.append(cmd)
            cmd = "chmod -R g+w {}".format(reference_fasta_fn + '.fai')
            cmds.append(cmd)
        cmd = 'longshot -P 0.01 -F -A --no_haps --bam {} --ref {} --out {} --potential_variants {}'.format(
            sorted_bam_fn, reference_fasta_fn, strand_bias_vcf_fn, merged_vcf_fn)
        cmds.append(cmd)
        cmd = "chmod -R g+w {}".format(strand_bias_vcf_fn)
        cmds.append(cmd)
    ## run annotation of merged vcf
    cmd = 'cd {} ; java -Xmx8g -jar snpEff.jar MN908947.3 {} >{}'.format(
        snpeff_dir, merged_vcf_fn, annotated_vcf_fn)
    cmds.append(cmd)
    cmd = "chmod -R g+w {}".format(annotated_vcf_fn)
    cmds.append(cmd)
    for cmd in cmds:
        retval = os.system(cmd)
        if retval != 0:
            logger.error('Command failed:\n{}'.format(cmd))
            exit(1)
    return True

def reformat_AA_change(s):
    # convert from three letter to single letter
    aa = "".join(aa_single_letter_.values()) + '*'
    for code, letter in aa_single_letter_.items():
        s = s.replace(code, letter)
    # convert AA changes
    m = re.fullmatch(f'p.([{aa}]+)(\d+)([{aa}]+)', s)
    if m:
        ref_aa, pos, alt_aa = m.group(1), int(m.group(2)), m.group(3)
        return ", ".join([f"{l}{pos+i}{r}" for i,(l,r) in enumerate(zip(ref_aa, alt_aa)) if l!=r])
    # convert deletions
    m = re.fullmatch(f'p.([{aa}]\d+_)?([{aa}]\d+)del', s)
    if m:
        if m.group(1) is not None:
            from_aa, to_aa = m.group(1)[:-1], m.group(2)
            if int(re.fullmatch(f"([{aa}])(\d+)", to_aa).group(2)) - int(re.fullmatch(f"([{aa}])(\d+)", from_aa).group(2)) > 1:
                return f"{from_aa}∆ - {to_aa}∆"
            else:
                return f"{from_aa}∆, {to_aa}∆"
        else:
            to_aa = m.group(2)
            return f"{to_aa}∆"
    # conservative inframe insertion
    m = re.fullmatch(f'p.([{aa}])(\d+)_([{aa}])(\d+)ins([{aa}]+)', s)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(1)}_{m.group(5)}"
        #return f"{m.group(3)}{m.group(4)}{m.group(5)}{m.group(3)}"
    # disruptive inframe insertion, p.D614delinsEL
    m = re.fullmatch(f'p.([{aa}])(\d+)delins([{aa}]+)', s)
    if m:
        return f"{m.group(1)}{m.group(2)}{m.group(3)}"
    # frame shift, p.P4715fs
    m = re.fullmatch(f'p.([{aa}])(\d+)fs', s)
    if m:
        return f"{m.group(1)}{m.group(2)}fs"
    return s

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

def is_masked(row, masked_regions):
    if pd.isnull(row[('final', 'decision')]) and not pd.isnull(row[('ARTIC', 'snv_filter')]):
        var_start = row[('medaka variant', 'site')]
        var_end = var_start + len(row[('medaka variant', 'ref')]) - 1
        if any((masked_regions['reference start'] <= var_start) & (masked_regions['reference end'] >= var_end)):
            row[('final', 'decision')] = 'masked'
        elif any((masked_regions['reference start'] <= var_start) & (masked_regions['reference end'] > var_start)) or \
             any((masked_regions['reference start'] < var_end) & (masked_regions['reference end'] >= var_end)) or \
             any((masked_regions['reference start'] > var_start) & (masked_regions['reference end'] <= var_end)):
            row[('final', 'decision')] = 'partially masked'
        else:
            row[('final', 'decision')] = 'rejected'
    return row

def enrich_introduced_variants(row):
    if pd.isnull(row[('medaka variant', 'site')]):
        row[('final', 'decision')] = 'introduced'
        m = re.fullmatch("(\D*)(\d+)(\D*)", str(row.name))
        if m is None:
            loger.error('variant id {} of unknown format'.format(row.name))
        ref, site, alt = m.group(1,2,3)
        row[('medaka variant', 'site')] = int(site)
        row[('medaka variant', 'ref')] = ref
        row[('medaka variant', 'alt')] = alt
    return row

def parse_ct_values(ct_values_fn):
    if os.path.exists(ct_values_fn):
        df = pd.read_csv(ct_values_fn, sep=';', header=0, names=['sample', 'ct'], dtype={'sample':str, 'ct':np.float32})
    else:
        logger.warning('File that should contain sample ct values does not exist: {}'.format(ct_values_fn))
        df = pd.DataFrame({'sample': pd.Series([], dtype='str'), 'ct':pd.Series([], dtype=np.float32)})
    return df.set_index('sample')

def load_snv_info(sample, artic_runs, results_dir, reference_fasta_fn, snpeff_dir, clades_df, subclades_df, alignment_tool):
    snv_info = []
    multiindex_rows = [[],[]]
    snpEff_fcts = {'Pool' : lambda info: info['Pool'],
                   "annotation": lambda info: info['ANN'][0].split("|")[1].replace("_", " ").replace("&", " & "),
                   "gene": lambda info: info['ANN'][0].split("|")[3],
                   "distance (nt)": lambda info: '-'+info['ANN'][0].split("|")[14] \
                   if info['ANN'][0].split("|")[1].startswith('upstream') \
                   else info['ANN'][0].split("|")[12].split("/")[0],
                   "impact": lambda info: info['ANN'][0].split("|")[2],
                   "AA change": lambda info: reformat_AA_change(info['ANN'][0].split("|")[10]),
                   "AA pos": lambda info: info['ANN'][0].split("|")[13]
                   }
    longshot_fcts = {
        'cov' : lambda info: int(info['AC'][0] + info['AC'][1] + info['AM']),
        '#ref' : lambda info: int(info['AC'][0]),
        '#alt' : lambda info: int(info['AC'][1]),
        '#amb' : lambda info: int(info['AM'])
    }
    artic_dir = artic_runs.loc[sample, 'artic_dir']
    annotated_vcf_fn = os.path.join(artic_dir, "{}.merged.snpeff.vcf".format(sample))
    pass_vcf_fn = os.path.join(artic_dir, "{}.pass.vcf.gz".format(sample))
    fail_vcf_fn = os.path.join(artic_dir, "{}.fail.vcf".format(sample))
    longshot_vcf_fn = os.path.join(artic_dir, "{}.longshot.vcf".format(sample))
    strand_bias_vcf_fn = os.path.join(artic_dir, "{}.longshot.01.vcf".format(sample))
    sample_final_consensus_fn = os.path.join(results_dir, sample, "{}.final.fasta".format(sample))
    # perform snv analysis if vcf files are missing
    if (not os.path.exists(annotated_vcf_fn)) or (not os.path.exists(strand_bias_vcf_fn)):
        if not run_extended_snv_pipeline(sample, artic_runs, snpeff_dir, reference_fasta_fn):
            logger.warning('Failed to run extended SNV analysis for sample {}'.format(sample))
    # parse vcf files
    vcf_ann = parse_vcf(annotated_vcf_fn, info_fcts=snpEff_fcts)
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
    vcf_annotation = vcf_ann[['annotation', 'gene', 'distance (nt)', 'impact', 'AA change', 'AA pos']]
    vcf_annotation.columns = pd.MultiIndex.from_product([['snpEff'], vcf_annotation.columns])
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
    # if a final consensus sequence is present in the results directory,
    # align it against the reference to extract information about SNVs that were confirmed or rejected
    if os.path.exists(sample_final_consensus_fn):
        vcf_confirmed, masked_regions, gap_start, gap_end = compare_consensus_to_reference(sample_final_consensus_fn, 
                                                                                           reference_fasta_fn, 
                                                                                           alignment_tool=alignment_tool)
        vcf_confirmed = vcf_confirmed[['decision', 'consensus site']]#.to_frame()
        vcf_confirmed.columns = pd.MultiIndex.from_product([['final'], vcf_confirmed.columns])
    else:
        vcf_confirmed = pd.DataFrame([], columns=pd.MultiIndex.from_product([['final'],['decision', 'consensus site']]))
        masked_regions, gap_start, gap_end = None, None, None

    # merge information from individual vcf files to one table
    if not vcf_medaka.loc[vcf_longshot.index.drop_duplicates()].index.equals(vcf_longshot.index):
        logger.warning('Number of variants of potentially different pools identified by ' + \
                       'medaka variant does not match the number of entries ' + \
                       'in the longshot output. Medaka variant information is therefore ' + \
                       'dropped for variants where an unambiguous assignment is not possible.')
    reconstruction = []
    for index in vcf_longshot.index.drop_duplicates():
        df_medaka = vcf_medaka.loc[index]
        if type(df_medaka) == pd.core.series.Series:
            df_medaka = df_medaka.to_frame().T
        df_longshot = vcf_longshot.loc[index]
        if type(df_longshot) == pd.core.series.Series:
            df_longshot = df_longshot.to_frame().T

        bias_count =  vcf_bias.index.to_list().count(index)
        if bias_count == len(df_longshot):
            df_longshot['strand bias'] = 'False'
        elif bias_count != 0:
            df_longshot['strand bias'] = 'Mixed'
        else:
            df_longshot['strand bias'] = 'True'

        df_longshot.columns = pd.MultiIndex.from_product([['longshot'], df_longshot.columns])
        if len(df_medaka) == len(df_longshot):
            reconstruction.append(pd.concat([df_medaka, df_longshot], axis=1))
        else:
            df_medaka = df_medaka.iloc[:len(df_longshot)]
            df_medaka[('medaka variant', 'Pool')] = 'N/A'
            df_medaka[('medaka variant', 'qual')] = np.nan
            reconstruction.append(pd.concat([df_medaka, df_longshot], axis=1))
    snv_info = pd.concat(reconstruction, axis=0, sort=False)

    columns = snv_info.columns #save column order, since pd.concat sorts the column names
    snv_info = pd.concat([snv_info, vcf_medaka.loc[vcf_medaka.index.difference(vcf_longshot.index.drop_duplicates())]]).sort_values(('medaka variant','site'))
    snv_info = snv_info.loc[:, columns] #restore column order
    snv_info = snv_info.join(vcf_artic_filter, how='left')
    snv_info = snv_info.join(vcf_annotation, how='left')
    snv_info = snv_info.join(vcf_confirmed, how='outer')

    # SNVs that were detected in the final fasta, but are not present in the previous list of potential SNVs
    sel = pd.isnull(snv_info[('medaka variant', 'site')])
    if np.any(sel):
        snv_info = snv_info.apply(enrich_introduced_variants, axis=1)

    # check if SNVs that are not present in the final fasta are masked with Ns
    if masked_regions is not None:
        snv_info = snv_info.apply(lambda row: is_masked(row, masked_regions), axis=1)

    # add nextstrain clade information
    clade_info = []
    for i,row in snv_info.iterrows():
        clades_ = []
        clades_.extend(list(clades_df.loc[(clades_df['site'] == row[('medaka variant', 'site')]) & \
                                          (clades_df['alt'] == row[('medaka variant', 'alt')]),'clade']))
        clades_.extend(list(subclades_df.loc[(subclades_df['site'] == row[('medaka variant', 'site')]) & \
                                             (subclades_df['alt'] == row[('medaka variant', 'alt')]),'clade']))
        clade_info.append(", ".join(clades_))
    snv_info.insert(loc=0, column=('nextstrain', 'clades'), value=clade_info)
    return snv_info.sort_values(('medaka variant', 'site')), masked_regions, gap_start, gap_end

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
    df = pd.Series(software_versions).to_frame()
    df.index.set_names(['software'], inplace=True)
    df.rename(columns={0: "version"}, inplace=True)
    return df
