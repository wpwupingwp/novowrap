#!/usr/bin/python3

from os import devnull
from pathlib import Path
from subprocess import run
from urllib.error import HTTPError
from urllib.request import urlopen
from zipfile import ZipFile
import argparse
import gzip
import logging

from Bio import SeqIO

from utils import get_fmt, get_ref, move
from validate_assembly import validate_main


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


def get_novoplasty():
    pl = list(Path('.').glob('NOVOPlasty*.pl'))
    if len(pl) != 0:
        return pl[0]
    _URL = 'https://github.com/ndierckx/NOVOPlasty/archive/NOVOPlasty3.6.zip'
    log.critical('Cannot find NOVOPlasty, try to download.')
    try:
        down = urlopen(_URL)
    except HTTPError:
        log.critical('Cannot download NOVOPlasty.')
        log.critical('Please manually download it from '
                     'http://github.com/ndierckx/NOVOPlasty')
        return None
    zip_file = Path('.') / 'NOVOPlasty3.6.zip'
    with open(zip_file, 'wb') as out:
        out.write(down.read())
    with ZipFile(zip_file, 'r') as z:
        # windows and linux both use "/"
        novoplasty = z.extract('NOVOPlasty-NOVOPlasty3.6/NOVOPlasty3.6.pl')
    zip_file.unlink()
    novoplasty = Path(novoplasty)
    log.info(f'Got {novoplasty.stem}.')
    return novoplasty


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    inputs = arg.add_argument_group('Input')
    inputs.add_argument('-f', help='forward fastq/gz file')
    inputs.add_argument('-r', help='reverse fastq/gz file')
    inputs.add_argument('-m', help='merged fastq/gz file')
    inputs.add_argument('-l', dest='list', help='csv file for batch mode')
    inputs.add_argument('-p', dest='platform', choices=['illumina', 'ion'],
                        default='illumina', help='sequencing platform')
    inputs.add_argument('-insert_size',
                        help='insert size of sequencing library')
    inputs.add_argument('-seed', default='rbcL,psaB,psaC,rrn23,rrn23s',
                        help='seed gene, separated by comma')
    inputs.add_argument('-seed_file',
                        help='seed file, will overwrite "-seed" option')
    inputs.add_argument('-split', default=0, type=int,
                        help='reads to use, set to 0 to skip split')
    options = arg.add_argument_group('Option')
    options.add_argument('-kmer', choices=range(23, 40, 2), default=39,
                         type=int, help='kmer size')
    options.add_argument('-min', default=100000, type=int,
                         help='minimum genome size (KB)')
    options.add_argument('-max', default=200000, type=int,
                         help='maximum genome size (KB)')
    options.add_argument('-mem', default=30, type=int,
                         help='maximum memory (GB)')
    options.add_argument('-o', dest='out', help='output folder')
    reference = arg.add_argument_group('Reference')
    reference.add_argument('-ref',
                           help='reference file, should be "gb" format with '
                           'only one record')
    reference.add_argument('-taxon', default='Nicotiana tabacum',
                           help='Taxonomy name')
    return arg.parse_args()


def check_arg(arg):
    """
    Check arg is ok or not.
    Set out if not give.
    Args:
        arg(NameSpace): arg generated by parse_args()
    Return:
        arg(NameSpace): arg with updated arguments
        arg_ok(bool) arg is valid or not
    """
    def _get_name(f, r):
        out_name = Path('Output').absolute()
        f = Path(f)
        r = Path(r)
        while f.suffix == r.suffix and f.suffix != '':
            f = f.with_suffix('')
            r = r.with_suffix('')
        same = 0
        idx = 0
        for i, j in zip(str(f), str(r)):
            if i == j:
                same += 1
            else:
                break
            idx += 1
        if same != 0:
            _ = list(str(f))
            _.pop(idx)
            strip_ = ''.join(_).rstrip('-_')
            if len(strip_) != 0:
                out_name = Path(strip_).absolute()
        return out_name

    if not any([arg.f, arg.r, arg.m, arg.list]):
        log.critical('Input is empty.')
        return arg, False
    if arg.out is None:
        if arg.f is not None and arg.r is not None:
            out = _get_name(arg.f, arg.r)
        elif arg.m is not None:
            out = Path(Path(arg.m).stem+'-out').absolute()
        elif arg.list is not None:
            out = Path(Path(arg.list).stem+'-out').absolute()
        else:
            out = Path('Output').absolute()
        arg.out = out
    else:
        arg.out = Path(arg.out).absolute()
    return arg, True


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
    """
    fmt = get_fmt(forward)
    new_f = output / Path(Path(forward).name).with_suffix(f'.{number}')
    new_r = output / Path(Path(reverse).name).with_suffix(f'.{number}')
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
    while count < number:
        # four line one record
        try:
            new_f_handle.write(next(f))
            new_f_handle.write(next(f))
            new_f_handle.write(next(f))
            new_f_handle.write(next(f))
            new_r_handle.write(next(r))
            new_r_handle.write(next(r))
            new_r_handle.write(next(r))
            new_r_handle.write(next(r))
        except StopIteration:
            break
        count += 1
    f_handle.close()
    r_handle.close()
    new_f_handle.close()
    new_r_handle.close()
    new_f = move(new_f, new_f.with_suffix(f'.{count}'))
    new_r = move(new_r, new_r.with_suffix(f'.{count}'))
    return new_f, new_r, count


def get_reads_len(filename):
    """
    Get reads length of fastq
    """
    fmt = get_fmt(filename)
    if fmt == 'gz':
        handle = gzip.open(filename)
    else:
        handle = open(filename, 'rb')
    handle.readline()
    seq = handle.readline()
    seq = seq.decode('utf-8').strip()
    length = len(seq)
    handle.close()
    return length


def get_seed(ref, output, gene):
    """
    Use BarcodeFinder to get seed or reference sequence.
    Arg:
        ref(Path): reference chloroplast genome gb file, only contains one
        record
        output(Path): output folder
        gene(str): gene names, separated by comma
    Return:
        seeds(list): seed files list
    """
    seeds = {}
    genes = gene.split(',')
    gb = SeqIO.read(ref, 'gb')
    accession = gb.annotations['accessions'][0]
    organism = gb.annotations['organism'][0].replace(' ', '_')
    for feature in gb.features:
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            gene_name = feature.qualifiers['gene'][0]
            # for rrn23S
            if gene_name in genes or gene_name.lower() in genes:
                seq = feature.extract(gb)
                seed_file = output / f'{gene_name}.fasta'
                with open(seed_file, 'w') as out:
                    out.write(f'>{gene_name}|{organism}|{accession}\n')
                    out.write(f'{seq.seq}\n')
                seeds[gene_name] = seed_file
    ordered_seeds = []
    for i in genes:
        if i in seeds:
            ordered_seeds.append(seeds[i])
    whole = ref.parent / 'whole.fasta'
    SeqIO.convert(ref, 'gb', whole, 'fasta')
    ordered_seeds.append(whole)
    return ordered_seeds


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
Read Length    = {arg.reads_len}
Insert size    = {arg.insert_size }
Platform       = {arg.platform}
Single/Paired  = {"PE" if arg.f is not None else "SE"}
Combined reads = {arg.m if arg.m is not None else ""}
Forward reads  = {arg.f if arg.f is not None else ""}
Reverse reads  = {arg.r if arg.r is not None else ""}

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
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
    # ensure novoplasty is available
    perl = run('perl -v', shell=True, stdout=open(devnull, 'w'))
    if perl.returncode != 0:
        log.critical('Please install Perl to run NOVOPlasty.')
        exit(-1)
    novoplasty = get_novoplasty()
    if novoplasty is None:
        exit(-1)
    arg = parse_args()
    arg, arg_ok = check_arg(arg)
    if not arg_ok:
        log.critical('Quit.')
        exit(-1)
    try:
        arg.out.mkdir()
    except FileExistsError:
        log.critical(f'Folder {arg.out.name} exists.')
        exit(-1)
    log_file_handler = logging.FileHandler(str(arg.out/'log.txt'))
    log_file_handler.setLevel(logging.INFO)
    Formatter = logging.Formatter(FMT, DATEFMT)
    log_file_handler.setFormatter(Formatter)
    log.addHandler(log_file_handler)
    log.info('Welcome to novowrap.')
    log.info(f'\tForward file: {arg.f}')
    log.info(f'\tReverse file: {arg.r}')
    log.info(f'\tSeeds: {arg.seed}')
    log.info(f'\tK-mer: {arg.kmer}')
    log.info(f'\tMinimum genome size: {arg.min}')
    log.info(f'\tMaximum genome size: {arg.max}')
    log.info(f'\tTaxonomy: {arg.taxon}')
    log.info(f'\tOutput folder: {arg.out}')
    if arg.f is not None:
        arg.reads_len = get_reads_len(arg.f)
    elif arg.m is not None:
        arg.reads_len = get_reads_len(arg.m)
    else:
        log.critical('Cannot detect reads length! Please input it manually.')
        arg.reads_len = int(input('Reads length:\t'))
    log.info(f'\tReads length: {arg.reads_len}')
    if arg.insert_size is None:
        arg.insert_size = arg.reads_len * 2 + 50
        log.info(f'The insert size is missing, use {arg.insert_size}.')
    if arg.split != 0:
        log.info(f'Split {arg.split} pairs of reads for assembly')
        arg.f, arg.r, splitted = split(arg.f, arg.r, arg.split, arg.out)
        if splitted < arg.split:
            log.warning(f'Want {arg.split} reads, acutally got {splitted}.')
    # get ref
    if arg.ref is not None:
        ref = Path(arg.ref)
        ref = move(ref, arg.out/ref, copy=True)
    else:
        log.info('Try to get reference from NCBI Genbank.')
        ref = get_ref(arg.taxon)
        if ref is None:
            log.critical('Cannot get reference.')
            exit(-1)
        else:
            log.info(f'Got {ref.stem}.')
            ref = move(ref, arg.out/ref)
    if arg.seed_file is not None:
        seeds = [Path(arg.seed_file), ]
    else:
        seeds = get_seed(ref, arg.out, arg.seed)
    if len(seeds) == 0:
        log.critical('Cannot get seeds!')
        exit(-1)
    csv_files = []
    success = False
    for seed in seeds:
        folder = arg.out / seed.stem
        folder.mkdir()
        log.info(f'Use {seed.stem} as seed.')
        config_file = config(arg.out, seed, arg)
        run(f'perl {novoplasty} -c {config_file}', shell=True)
        # novoplasty use current folder as output folder
        circularized, options, merged, contigs = organize_out(
            Path().cwd(), folder)
        if len(circularized) == 0 and len(options) == 0 and len(merged) == 0:
            log.warning(f'Assembled with {seed.stem} failed.')
            continue
        validated = []
        # validate merged or not?
        log.info('Validate assembly results.')
        # for i in (*circularized, *options, *merged):
        for i in (*circularized, *options):
            arg_str = f'{i} -ref {ref} -seed {seed.stem} -o {folder}'
            validate_file, report = validate_main(arg_str)
            validated.extend(validate_file)
            if report not in csv_files:
                csv_files.append(report)
        if len(validated) != 0:
            success = True
            break
        else:
            log.warning('No records passed validation.')
        if not success:
            log.critical(f'Assembly with {seed.stem} failed.')
    if len(csv_files) != 0:
        merged_csv = move(csv_files[0], arg.out / 'Validation.csv', copy=True)
        if len(csv_files) > 1:
            with open(merged_csv, 'a') as h1:
                for i in csv_files[1:]:
                    with open(i, 'r') as h2:
                        h1.write(''.join(h2.readlines()[1:]))
    log.info(f'Merged validation results were written into {merged_csv.name}')
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
