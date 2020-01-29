#!/usr/bin/python3

from os import chdir
from pathlib import Path
from random import randint
from subprocess import DEVNULL, run
from threading import Thread
from time import sleep
from urllib.request import urlopen
from zipfile import ZipFile
import argparse
import gzip
import logging

from Bio import SeqIO

from utils import get_fmt, get_ref, accessible, move
from merge import merge_main
from validate import validate_main


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


def get_novoplasty(arg):
    """
    Ensure perl and novoplasty is available.
    Return novoplasty's path or None.
    Args:
        arg(NameSpace): args
    Return:
        novoplasty(Path): path of perl file
    """
    url = 'https://github.com/ndierckx/NOVOPlasty/archive/NOVOPlasty3.7.2.zip'
    perl = run('perl -v', shell=True, stdout=DEVNULL, stderr=DEVNULL)
    if perl.returncode != 0:
        log.critical('Please install Perl for running NOVOPlasty.')
        return None
    filename = 'NOVOPlasty3.7.2.pl'
    pl = Path('.').absolute() / filename
    if pl.exists():
        log.debug('Found NOVOPlasty in current folder.')
        return pl
    pl = arg.third_party / filename
    if pl.exists():
        log.debug('Found NOVOPlasty in third_party folder.')
        return pl
    log.critical('Cannot find NOVOPlasty, try to download.')
    log.info('\tThe program assumes that users accept the license of '
             'NOVOPlasty. ')
    log.info('\tSee https://raw.githubusercontent.com/'
             'ndierckx/NOVOPlasty/ for details.')
    try:
        # 126kb/10s=12.6kb/s, enough for test
        _ = urlopen('https://github.com', timeout=10)
    except Exception:
        log.critical('Cannot connect to github.com.')
        log.critical('Please check your Internet connection.')
        return None
    log.info('Due to connection speed, may need minutes. Please be patient.')
    try:
        # file is ~18mb, 18mb/360s=50kb/s, consider it's ok for most of users
        down = urlopen(url, timeout=360)
    except Exception:
        log.critical('Cannot download NOVOPlasty.')
        log.critical('Please manually download it from '
                     'https://github.com/ndierckx/NOVOPlasty')
        return None
    zip_file = arg.third_party / 'NOVOPlasty3.7.2.zip'
    with open(zip_file, 'wb') as out:
        out.write(down.read())
    with ZipFile(zip_file, 'r') as z:
        # windows and linux both use "/"
        novoplasty = z.extract('NOVOPlasty-NOVOPlasty3.7.2/NOVOPlasty3.7.2.pl')
    novoplasty = Path(novoplasty).absolute()
    novoplasty = move(novoplasty, pl)
    log.info(f'Got {novoplasty.stem}.')
    return novoplasty


def parse_args(arg_list=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    inputs = arg.add_argument_group('Input')
    inputs.add_argument('-i', dest='input', nargs='*',
                        help='single or pair-end input, fastq or gz format')
    inputs.add_argument('-l', dest='list', help='csv file for batch mode')
    inputs.add_argument('-p', dest='platform', choices=['illumina', 'ion'],
                        default='illumina', help='sequencing platform')
    inputs.add_argument('-insert_size',
                        help='insert size of sequencing library')
    inputs.add_argument('-seed', default='rbcL,psaB,psaC,rrn23',
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
    options.add_argument('-o', '-out', dest='out', help='output folder')
    options.add_argument('-debug', action='store_true', help='debug mode')
    reference = arg.add_argument_group('Reference')
    reference.add_argument('-ref',
                           help='reference file, should be "gb" format with '
                           'only one record')
    reference.add_argument('-taxon', default='Nicotiana tabacum',
                           help='Taxonomy name')
    validate = arg.add_argument_group('Validate')
    validate.add_argument('-perc_identity', type=float, default=0.7,
                          help='minimum percentage of identity of BLAST, 0-1')
    validate.add_argument('-len_diff', type=float, default=0.2,
                          help='maximum percentage of length differnce of '
                          'query to reference, 0-1')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def get_output(arg):
    """
    Get output folder.
    If exists, return None.
    Args:
        arg(NameSpace): arg generated by parse_args()
    Return:
        out(Path or None): output path
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

    out = Path('.')
    if arg.out is None:
        if arg.input is not None:
            if len(arg.input) == 2:
                out = _get_name(arg.input[0], arg.input[1])
            else:
                # for single file, directly remove all suffixes is dangerous
                out = Path(Path(arg.input[0]).stem).absolute()
        elif arg.list is not None:
            out = Path(Path(arg.list).stem).absolute()
    else:
        out = Path(arg.out).absolute()
    if out.exists():
        log.warning(f'Output folder {out.name} exists.')
        # give users one more chance
        new_name = out.name + f'-{randint(0, 2020):04d}'
        out = out.with_name(new_name)
        if out.exists():
            log.critical('Cannot create output folder.')
            return None
        else:
            log.info(f'Use {out.name} instead.')
    return out.absolute()


def init_arg(arg):
    """
    Initialize working folder with arg.
    Args:
        arg(NameSpace): arg generated from parse_args
    Return:
        success(bool): success or not
        arg(NameSpace): arg with output info
    """
    if arg.debug:
        logging.basicConfig(level=logging.DEBUG)
        try:
            import coloredlogs
            coloredlogs.install(level=logging.DEBUG, fmt=FMT, datefmt=DATEFMT)
        except ImportError:
            pass
    success = False
    if arg.list is None and arg.input is None:
        log.critical('Input is empty.')
        return success, arg
    if len(arg.input) > 2:
        log.critical('Only accept one or two input file(s).')
        return success, arg
    if arg.list is not None:
        if not arg.list.exists():
            log.critical(f'Input file {arg.list} does not exists.')
            return success, arg
        arg.list = Path(arg.list).absolute()
    if arg.input is not None:
        arg.input = [Path(i).absolute() for i in arg.input]
        for i in arg.input:
            if not i.exists():
                log.critical(f'Input file {i} does not exists.')
                return success, arg
    arg.out = get_output(arg)
    if arg.out is None:
        return success, arg
    if not accessible(arg.out, 'folder'):
        log.critical(f'Failed create {arg.out}. Please contact the '
                     f'administrator.')
        return success, arg
    arg.out.mkdir()
    arg.log = arg.out / 'Log'
    arg.log.mkdir()
    arg.raw = arg.out / 'Raw'
    arg.raw.mkdir()
    arg.tmp = arg.out / 'Temp'
    arg.tmp.mkdir()
    arg.third_party = Path().home().absolute() / '.novowrap'
    if not arg.third_party.exists():
        log.debug(f'Create folder {arg.third_party}')
        arg.third_party.mkdir()
    if not accessible(arg.third_party/'test', 'file'):
        log.critical(f'Failed to access {arg.third_party}.'
                     f'Please contact the administrator.')
        return success, arg
    success = True
    return success, arg


def read_table(arg):
    """
    Read table from given csv file.
    Columns of table:
        Input, Input(optional), Taxonomy
    Args:
        arg(NameSpace): arg generated from parse_args
    Return:
        inputs(list): [[f, r], taxon]
    """
    inputs = []
    with open(arg.list, 'r', encoding='utf-8') as raw:
        for line in raw:
            try:
                f, r, taxon = line.strip().split(',')
            except IndexError:
                log.warning(f'Cannot parse the line : {line}')
                continue
            if r == '':
                f_r = [f, ]
            else:
                f_r = [f, r]
            inputs.append([f_r, taxon])
    return inputs


def split(raw, number, output):
    """
    Split reads of original file from the beginning.
    Args:
        raw(str or Path): input file, coulde be fastq or gz format
        number(int): number of reads to split
        output(Path): output folder
    Return:
        splitted(Path): splitted file, fastq format
        count(int): reads actually got
    """
    raw = Path(raw).absolute()
    fmt = get_fmt(raw)
    splitted = output / raw.with_suffix(f'.{number}').name
    splitted_handle = open(splitted, 'wb')
    if fmt == 'gz':
        raw_handle = gzip.open(raw)
    else:
        raw_handle = open(raw, 'rb')
    line = iter(raw_handle)
    count = 0
    while count < number:
        # four line one record
        try:
            splitted_handle.write(next(line))
            splitted_handle.write(next(line))
            splitted_handle.write(next(line))
            splitted_handle.write(next(line))
        except StopIteration:
            break
        count += 1
    raw_handle.close()
    splitted_handle.close()
    splitted = move(splitted, splitted.with_suffix(f'.{count}'))
    return splitted, count


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
    log.debug(f'\tReads length: {length}')
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
    organism = gb.annotations['organism'].replace(' ', '_')
    for feature in gb.features:
        if feature.type == 'gene' and 'gene' in feature.qualifiers:
            gene_name = feature.qualifiers['gene'][0]
            if gene_name in genes:
                seq = feature.extract(gb)
                seed_file = output / f'{gene_name}.seed'
                with open(seed_file, 'w', encoding='utf-8') as out:
                    out.write(f'>{gene_name}|{organism}|{accession}\n')
                    out.write(f'{seq.seq}\n')
                seeds[gene_name] = seed_file
    ordered_seeds = []
    for i in genes:
        if i in seeds:
            ordered_seeds.append(seeds[i])
    whole = output / 'whole.seed'
    SeqIO.convert(ref, 'gb', whole, 'fasta')
    ordered_seeds.append(whole)
    return ordered_seeds


def config(seed, arg):
    """
    Generate config file for NOVOPlasty.
    Arg:
        seed(Path): seed file
        arg(NameSpace): parameters user provided
    Return:
        config_file(Path): config file
    """
    if len(arg.input) == 2:
        f, r = arg.input
        m = ''
        s_or_p = 'PE'
    else:
        f = r = ''
        m = arg.input[0]
        s_or_p = 'SE'
    arg.reads_len = get_reads_len(arg.input[0])
    if arg.insert_size is None:
        arg.insert_size = arg.reads_len * 2 + 50
        log.debug(f'The insert size is missing, use {arg.insert_size}.')
    config_str = f"""Project:
-----------------------
Project name          = {arg.out.name}
Type                  = chloro
Genome Range          = {arg.min}-{arg.max}
K-mer                 = {arg.kmer}
Max memory            = {arg.mem}
Extended log          = 1
Extend seed directly  = no
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
Single/Paired  = {s_or_p}
Combined reads = {m}
Forward reads  = {f}
Reverse reads  = {r}

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.9
Insert Range strict   = 1.3
Use Quality Scores    = no
"""
    config_file = arg.raw / f'{seed.stem}_config.ini'
    with open(config_file, 'w', encoding='utf-8') as out:
        out.write(config_str)
    return config_file


def organize_out(arg, seed):
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
        arg(NameSpace): args
        seed(str): seed gene's name
    Return:
        contigs(list): contig files
        merged(list): merged files
        options(list): options files
        circularized(list): circularized files
    """
    def txt_to_fasta(old):
        """
        Convert NOVOPlasty-generated txt to standard fasta.
        """
        clean = []
        record = []
        begin = False
        with open(old, 'r', encoding='utf-8') as raw:
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
        with open(new, 'w', encoding='utf-8') as output:
            for line in clean:
                output.write(line.replace('*', ''))
        return new

    for i in arg.out.glob('contigs_tmp_*.txt'):
        i = i.absolute()
        move(i, arg.tmp/i.with_name(f'{i.stem}-{seed}{i.suffix}').name)
    for i in arg.out.glob('log_*.txt'):
        i = i.absolute()
        move(i, arg.log/i.with_name('NOVOPlasty-'+i.name).name)
    contigs = []
    for i in arg.out.glob('Contigs_*.fasta'):
        i = i.absolute()
        i = move(i, arg.raw/i.with_name(f'{i.stem}-{seed}{i.suffix}').name)
        contigs.append(i)
    merged = []
    for i in arg.out.glob('Merged_contigs_*.txt'):
        i = i.absolute()
        i = move(i, arg.raw/i.with_name(f'{i.stem}-{seed}{i.suffix}').name)
        fasta = txt_to_fasta(i)
        fasta = move(fasta, arg.raw/fasta.name)
        merged.append(fasta)
    options = []
    for i in arg.out.glob('Option_*.fasta'):
        i = i.absolute()
        i = move(i, arg.raw/i.with_name(f'{i.stem}-{seed}{i.suffix}').name)
        options.append(i)
    circularized = []
    for i in arg.out.glob('Circularized_assembly*.fasta'):
        i = i.absolute()
        i = move(i, arg.raw/i.with_name(f'{i.stem}-{seed}{i.suffix}').name)
        circularized.append(i)
    return circularized, options, merged, contigs


def assembly(arg, novoplasty):
    """
    Assembly input file by wrapping NOVOPlasty.
    Args:
        arg(NameSpace): arguments
        novoplasty(Path): novoplasty file
    Return:
        success(Bool): success or not
    """
    def _patient_log():
        # hint user every 30s to avoid long time boring waiting
        n = 0
        while novoplasty_is_running:
            sleep(1)
            n += 1
            if n >= 30:
                log.info('NOVOPlasty is running, please be patient...')
                n = 0
        return

    success = False
    log.info('')
    for i in arg.input:
        test = Path(i)
        # arg.list may contains invalid file
        if not test.exists():
            log.critical(f'Cannot find input file {i}')
            return success
        else:
            log.info(f'Input file: {i}')
    log.info(f'Minimum genome size: {arg.min}')
    log.info(f'Maximum genome size: {arg.max}')
    log.info(f'Taxonomy: {arg.taxon}')
    log.info(f'Output folder: {arg.out}')
    if arg.split != 0:
        log.info(f'Split {arg.split} pairs of reads for assembly')
        splitted = []
        for raw in arg.input:
            new, count = split(raw, arg.split, arg.tmp)
            if count < arg.split:
                log.warning(f'Want {arg.split} reads, acutally got {count}.')
            splitted.append(new)
        arg.input = splitted
    # get ref
    if arg.ref is not None:
        if get_fmt(arg.ref) != 'gb':
            log.critical('Reference file should be genbank format, '
                         'but {arg.ref} is not.')
            return success
        ref = Path(arg.ref).absolute()
        ref = move(ref, arg.tmp/ref.name, copy=True)
    else:
        ref, arg.taxon = get_ref(arg.taxon, arg.tmp)
        if ref is None:
            log.critical('Cannot get reference.')
            return success
        else:
            ref = move(ref, arg.tmp/ref.name)
    # get seed
    seeds = []
    ordered_seeds = get_seed(ref, arg.raw, arg.seed)
    if arg.seed_file is not None:
        seeds.append(Path(arg.seed_file).absolute())
        # only add whole.seed
        seeds.append(ordered_seeds[-1])
    else:
        seeds.extend(ordered_seeds)
    if len(seeds) == 0:
        log.critical('Cannot get seeds!')
        return success
    csv_files = []
    all_contigs = []
    for seed in seeds:
        log.info(f'Use {seed.stem} as seed.')
        config_file = config(seed, arg)
        log.info('Call NOVOPlasty... May need minutes (rarely half an hour)')
        # use mark to terminate thread
        novoplasty_is_running = True
        hint = Thread(target=_patient_log)
        hint.start()
        # ignore bad returncode
        run(f'perl {novoplasty} -c {config_file}', shell=True,
            stdout=DEVNULL, stderr=DEVNULL)
        novoplasty_is_running = False

        # novoplasty use current folder as output folder
        circularized, options, merged, contigs = organize_out(arg, seed.stem)
        all_contigs.extend(contigs)
        if len(circularized) == 0 and len(options) == 0 and len(merged) == 0:
            log.warning(f'Assembled with {seed.stem} failed.')
            continue
        validated = []
        log.info('Validate assembly results.')
        # validate merged or not?
        for i in (*circularized, *options):
            arg_str = (f'{i} -ref {ref} -seed {seed.stem} -o {arg.out} '
                       f'-perc_identity {arg.perc_identity} '
                       f'-len_diff {arg.len_diff}')
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
    if not success:
        log.info('Failed with all seeds.')
        all_input = ' '.join([str(i.absolute()) for i in all_contigs])
        # '' means empty
        if len(all_input) == 0:
            return success
        log.info('Try to assembly contigs generated from each seed.')
        arg_str = f'{all_input} -o {arg.out/"Raw"/"merge_seed.fasta"}'
        n_assembly, assembly_result = merge_main(arg_str)
        if n_assembly != 0:
            arg_str = (f'{assembly_result} -ref {ref} -seed merge '
                       f'-o {arg.out} '
                       f'-perc_identity {arg.perc_identity} '
                       f'-len_diff {arg.len_diff}')
            validate_file, report = validate_main(arg_str)
            if len(validate_file) != 0:
                success = True
                csv_files.append(report)
    return success


def assembly_main(arg_str=None):
    """
    For consistence, return (success, arg.out).
    """
    log.info('Welcome to novowrap.')
    # check arg
    if arg_str is None:
        arg = parse_args()
    else:
        arg = parse_args(arg_str.split(' '))
    success, arg = init_arg(arg)
    if not success:
        log.critical('Quit.')
        return success, arg.out
    else:
        log.debug('Init OK.')
    cwd = Path().cwd().absolute()
    chdir(arg.out)
    # check before run
    novoplasty = get_novoplasty(arg)
    if novoplasty is None:
        log.critical('Quit.')
        success = False
        return success, arg.out
    # log to file
    log_file_handler = logging.FileHandler(str(arg.log/'Log.txt'))
    # more detail in file log
    log_file_handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter(FMT, DATEFMT)
    log_file_handler.setFormatter(formatter)
    log.addHandler(log_file_handler)
    # start
    if arg.list is None:
        success = assembly(arg, novoplasty)
    else:
        table = read_table(arg)
        for i in table:
            arg.input, arg.taxon = i
            success = assembly(arg, novoplasty)
    log.info('Bye.')
    chdir(cwd)
    return success, arg.out


if __name__ == '__main__':
    assembly_main()
