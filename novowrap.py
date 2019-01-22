#!/usr/bin/python3

from pathlib import Path
from subprocess import run
from timeit import default_timer as timer
import argparse


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('-f', required=True, help='forward fastq file')
    arg.add_argument('-r', required=True, help='reverse fastq file')
    arg.add_argument('-kmer', choices=range(23, 40, 2), default=39, type=int,
                     help='kmer size')
    arg.add_argument('-min', default=120_000, help='minimum genome size (kb)')
    arg.add_argument('-max', default=200_000, help='minimum genome size (kb)')
    arg.add_argument('-reads_len', default=150, help='reads length')
    arg.add_argument('-seed', help='seed file')
    # arg.add_argument('-split', default=1_000_000,
    #                  help='reads to use (million), set to 0 to skip split')
    return arg.parse_args()


def get_seed(arg):
    folder = Path('seed')
    seed_candidate = [folder/i for i in ('rbcL.fasta', 'rrn23.fasta',
                                         'psaC.fasta', 'rrn5.fasta')]
    if arg.seed is not None:
        seed_candidate = [arg.seed, *seed_candidate]
    for i in seed_candidate:
        if i not in arg.seed_failed:
            return i
    raise Exception('All seeds failed.')


def config(name, arg):
    template = """Project:
-----------------------
Project name          = {name}
Type                  = chloro
Genome Range          = {min_size}-{max_size}
K-mer                 = {kmer}
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = {seed}
Reference sequence    =
Variance detection    = no
Heteroplasmy          =
HP exclude list       =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = {reads_len}
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = {forward}
Reverse reads         = {reverse}

Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.8
Insert Range strict   = 1.3
Use Quality Scores    = no
"""
    config = template.format(name=name, min_size=arg.min, max_size=arg.max,
                             kmer=arg.kmer, seed=arg.seed,
                             reads_len=arg.reads_len, forward=arg.f,
                             reverse=arg.r)
    config_file = name / 'config.ini'
    with open(config_file, 'w') as out:
        out.write(config)
    return config_file


def clean(name):
    fq = name.stem
    contigs = Path('.').glob('Contigs_*{}*'.format(fq))
    options = Path('.').glob('Option_*{}*'.format(fq))
    merged = Path('.').glob('Merged_contigs_{}.txt'.format(fq))
    tmp = Path('.').glob('contigs_tmp_{}.txt'.format(fq))
    log = Path('.').glob('log_{}.txt'.format(fq))
    for i in [*contigs, *options, *merged, *tmp, *log]:
        i.rename(name/i)


def main():
    start = timer()
    arg = parse_args()
    arg.seed_failed = {}
    name = Path(Path(arg.f).stem)
    name.mkdir()
    while True:
        arg.seed = get_seed(arg)
        config_file = config(name, arg)
        test = run('perl NOVOPlasty2.7.2.pl -c {}'.format(config_file),
                   shell=True)
        if test.returncode == 0:
            clean(name)
            break
        else:
            arg.seed_failed.add(arg.seed)
        return -1
    end = timer()
    print('Cost {:.3f} secondes.'.format(end-start))
    return


if __name__ == '__main__':
    main()
