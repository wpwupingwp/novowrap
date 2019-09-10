#!/usr/bin/python3

from os import devnull
from pathlib import Path
from subprocess import run
import argparse
import gzip
import logging

from Bio import SeqIO

from utils import get_fmt, get_ref, move
from validate_assembly import validate_main


NULL = open(devnull, 'w')
# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('novowrap')
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    inputs = arg.add_argument_group('Input')
    inputs.add_argument('-f', required=True, help='forward fastq/gz file')
    inputs.add_argument('-r', required=True, help='reverse fastq/gz file')
    inputs.add_argument('-split', default=0, help='reads to use, set to 0 '
                        'to skip split')
    options = arg.add_argument_group('Option')
    options.add_argument('-kmer', choices=range(23, 40, 2), default=39,
                         type=int, help='kmer size')
    options.add_argument('-min', default=100000, type=int,
                         help='minimum genome size (KB)')
    options.add_argument('-max', default=200000, type=int,
                         help='maximum genome size (KB)')
    options.add_argument('-mem', default=30, type=int,
                         help='maximum memory (GB)')
    options.add_argument('-reads_len', default=150, type=int,
                         help='reads length')
    options.add_argument('-gene', help='seed gene')
    options.add_argument('-try', dest='try_n', type=int,
                         default=5, help='maximum tried times')
    ref = arg.add_argument_group('Reference')
    ref.add_argument('-taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
    return arg.parse_args()


def split(forward, reverse, number, output):
    """
    Split reads of original file from the beginning.
    Args:
        forward(str or Path): forward file, could be fastq or gz
        reverse(str or Path): reverse file, could be fastq or gz
        number(int): number of reads to split
        output(Path): output folder
    Return:
        new_f(Path): new forward file
        new_r(Path): new reverse file
        count(n): number of reads got, may smaller than target number
    """
    fmt = get_fmt(forward)
    new_f = output / Path(forward).with_suffix(f'.{number}')
    new_r = output / Path(reverse).with_suffix(f'.{number}')
    new_f_handle = open(new_f, 'wb')
    new_r_handle = open(new_r, 'wb')
    if fmt == 'gz':
        f_handle = gzip.open(forward)
        r_handle = gzip.open(reverse)
    else:
        f_handle = open(forward, 'rb')
        r_handle = open(reverse, 'rb')
    f = iter(f_handle)
    r = iter(r_handle)
    count = 0
    while count <= number:
        # four line one record
        new_f_handle.write(next(f))
        new_f_handle.write(next(f))
        new_f_handle.write(next(f))
        new_f_handle.write(next(f))
        new_r_handle.write(next(r))
        new_r_handle.write(next(r))
        new_r_handle.write(next(r))
        new_r_handle.write(next(r))
        count += 1
    f_handle.close()
    r_handle.close()
    new_f_handle.close()
    new_r_handle.close()
    return new_f, new_r, count


def get_seed(ref, output, gene=None):
    """
    Use BarcodeFinder to get seed or reference sequence.
    Arg:
        ref(Path): reference chloroplast genome gb file, only contains one
        record
        output(Path): output folder
        gene(tuple): gene name
    Return:
        seeds(list): seed files list
    """
    # strand: +, -, -, -, +
    candidate_genes = ['rbcL', 'matK', 'psaB', 'psaC', 'rrn23']
    seeds = []
    if gene is None:
        genes = candidate_genes
    else:
        genes = [gene, ]
        genes.extend(candidate_genes)
    gb = SeqIO.read(ref, 'gb')
    accession = gb.annotations['accessions'][0]
    organism = gb.annotations['organism'][0].replace(' ', '_')
    for feature in gb.features:
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            gene_name = feature.qualifiers['gene'][0]
            if gene_name in genes:
                seq = feature.extract(gb)
                file = output / f'{gene_name}.fasta'
                with open(file, 'w') as out:
                    out.write(f'>{gene_name}|{organism}|{accession}\n')
                    out.write(f'{seq}\n')
                seeds.append(file)
    return seeds


def config(out, seed, arg):
    """
    Generate config file for NOVOPlasty.
    Arg:
        out(Path): output folder
        seed(Path): seed file
        arg(NameSpace): parameters user provided
    Return:
        config_file(Path): config file
    """
    config = f"""Project:
-----------------------
Project name          = {out.name}
Type                  = chloro
Genome Range          = {arg.min}-{arg.max}
K-mer                 = {arg.kmer}
Max memory            = {arg.mem}
Extended log          = 1
Save assembled reads  = no
Seed Input            = {seed}
Reference sequence    =
Variance detection    = no
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = {arg.reads_len}
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = {arg.f}
Reverse reads         = {arg.r}

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no
"""
    config_file = out / f'{seed.stem}_config.ini'
    with open(config_file, 'w') as out:
        out.write(config)
    return config_file


def txt_to_fasta(old):
    clean = []
    record = []
    begin = False
    with open(old, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                clean.extend(record)
                record = []
                begin = True
            if line.startswith(' ') or len(line.strip()) == 0:
                begin = False
                clean.extend(record)
                record = []
            if begin:
                record.append(line)
    clean.extend(record)
    new = Path(old).with_suffix('.fasta')
    with open(new, 'w') as out:
        for line in clean:
            out.write(line.replace('*', ''))
    return new


def organize_out(source, dest):
    """
    Organize NOVOPlasty output.
        log*: log file
        contigs_tmp*: temporary files
        Contigs*: contigs
        Merged*: merged contigs, may be circular or empty, contains options
        Option*: merged contigs, circular or incomplete circular
        Circularized*: circularized sequence
    Return fasta list.
    Arg:
        source(Path): current directory
        dest(Path): directory to move
    Return:
        contigs(list): contig files
        merged(list): merged files
        options(list): options files
        circularized(list): circularized files
    """
    for i in source.glob('contigs_tmp_*'):
        i = move(i, dest/i.name)
    for i in source.glob('log_*.txt'):
        i = move(i, dest/i.name)
    contigs = [txt_to_fasta(move(i, dest/i.name)) for i in
               source.glob('Contigs_*')]
    merged = [txt_to_fasta(move(i, dest/i.name)) for i in
              source.glob('Merged_contigs_*')]
    options = [txt_to_fasta(move(i, dest/i.name)) for i in
               source.glob('Option_*')]
    circularized = [txt_to_fasta(move(i, dest/i.name)) for i in
                    source.glob('Circularized_assembly*')]
    return circularized, options, merged, contigs


def main():
    arg = parse_args()
    out = Path(Path(arg.f).stem+'-out').absolute()
    try:
        out.mkdir()
    except FileExistsError:
        log.critical(f'Folder {out.name} exists.')
        exit(-1)
    log_file_handler = logging.FileHandler(str(out/'log.txt'))
    log_file_handler.setLevel(logging.INFO)
    Formatter = logging.Formatter(FMT, DATEFMT)
    log_file_handler.setFormatter(Formatter)
    log.addHandler(log_file_handler)
    log.info('Welcome to novowrap.')
    log.info(f'Forward file:\t{arg.f}')
    log.info(f'Reverse file:\t{arg.r}')
    log.info(f'K-mer:\t{arg.kmer}')
    log.info(f'Minimum genome size:\t{arg.min}')
    log.info(f'Maximum genome size:\t{arg.max}')
    log.info(f'Reads length:\t{arg.reads_len}')
    log.info(f'Taxonomy:\t{arg.taxon}')
    log.info(f'Maximum tried times:\t{arg.try_n}')
    log.info(f'Use {out} as output folder.')
    if arg.split != 0:
        log.info(f'Split {arg.split} pairs of reads for assembly')
        arg.f, arg.r = split(arg.f, arg.r, arg.split, out)
    success = False
    fail = 0
    taxon, ref, accession = get_ref(arg.taxon)
    if taxon is None:
        exit(-1)
    ref = move(ref, out/ref)
    seeds = get_seed(ref, out, arg.gene)
    if len(seeds) == 0:
        log.critical('Cannot get seeds!')
        exit(-1)
    for seed in seeds:
        if fail >= arg.try_n:
            log.critical(f'Too much failure. Quit.')
            break
        folder = out / seed.stem
        folder.mkdir()
        log.info(f'No. {fail+1} try, use {seed.name} as seed file.')
        config_file = config(out, seed, arg)
        log.info('NOVOPlasty version:\t3.6')
        run_novo = run(f'perl NOVOPlasty3.6.pl -c {config_file}', shell=True)
        if run_novo.returncode != 0:
            log.critical('Failed to run NOVOPlasty. Quit.')
            exit(-1)
        # novoplasty generates outputs in current folder
        # use rbcL to detect strand direction
        log.info(f'Organize NOVOPlasty output of {seed.name}.')
        # novoplasty use current folder as output folder
        circularized, options, merged, contigs = organize_out(
            Path().cwd(), folder)
        if len(circularized) == 0 and len(options) == 0 and len(merged) == 0:
            log.warning(f'Assembled with {seed.name} failed.')
            fail += 1
            continue
        validated = []
        # validate merged or not?
        log.info('Validate assembly results.')
        for i in (*circularized, *options, *merged):
            arg_str = f'{i} -r {ref} -o {folder}'
            validated.append(validate_main(arg_str))
        if len(validated) != 0:
            success = True
            break
        else:
            log.warning('Cannot find correct conformation for all assembly of '
                        f'{seed.name}.')
            fail += 1
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    NULL.close()
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
