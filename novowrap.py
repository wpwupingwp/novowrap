#!/usr/bin/python3

from Bio import Entrez
from pathlib import Path
from platform import system
from subprocess import run
import argparse
import logging
import re


# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%I:%M:%S'
TEMP_LOG = 'Temp.log'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO,
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler(TEMP_LOG)])
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
log = logging.getLogger(__name__)


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('-f', required=True, help='forward fastq file')
    arg.add_argument('-r', required=True, help='reverse fastq file')
    arg.add_argument('-kmer', choices=range(23, 40, 2), default=39, type=int,
                     help='kmer size')
    arg.add_argument('-min', default=40000, help='minimum genome size (kb)')
    arg.add_argument('-max', default=300000, help='minimum genome size (kb)')
    arg.add_argument('-reads_len', default=150, help='reads length')
    arg.add_argument('-taxon', required=True, help='Taxonomy name')
    # arg.add_argument('-split', default=1_000_000,
    #                  help='reads to use (million), set to 0 to skip split')
    return arg.parse_args()


def get_full_taxon(taxon):
    """
    Get string.
    Return lineage list, only contains Kingdom, Phylum, Class, Order, Family,
    Genus, Species.
    """
    split = taxon.split(' ')
    if len(split) >= 2 and split[0][0].isupper() and split[0][1].islower():
        name = ' '.join(split[0:2])
    else:
        name = split[0]
    Entrez.email = 'guest@example.org'
    search = Entrez.read(Entrez.esearch(db='taxonomy', term=f'"{name}"'))
    if search['Count'] == '0':
        log.critical(f'Cannot find {name} in NCBI Taxonomy.')
    taxon_id = search['IdList'][0]
    record = Entrez.read(Entrez.efetch(db='taxonomy', id=taxon_id))[0]
    names = [i['ScientificName'] for i in record['LineageEx']]
    full_lineage = {i['Rank']: i['ScientificName'] for i in
                    record['LineageEx']}
    full_lineage[record['Rank']] = record['ScientificName']
    if 'kingdom' not in full_lineage:
        full_lineage['kingdom'] = full_lineage['superkingdom']
    if ('class' not in full_lineage and 'order' in full_lineage):
        last_phyta = ''
        for i in names[::-1]:
            if i.endswith('phyta'):
                last_phyta = i
                break
        # virus do not have phylum?
        phylum = full_lineage.get('phylum', None)
        if last_phyta != phylum:
            full_lineage['class'] = last_phyta
    target = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus',
              'species']
    lineage = []
    for i in target:
        lineage.append(full_lineage.get(i, ''))
    if ' ' in lineage[-1]:
        lineage[-1] = lineage[-1].split(' ')[-1]
    return lineage


def get_seed(taxon, output):
    genes = ('rbcL', 'matK', 'psaB', 'psaC', 'rrn23', 'rrn5')
    MAX_LEN = 250000
    if system() == 'Windows':
        python = 'python'
    else:
        python = 'python3'
    for taxon in reversed(get_full_taxon(taxon)):
        if taxon == '':
            continue
        for gene in genes:
            out = output / f'{taxon}-{gene}'
            log.info(f'Querying {out}.')
            # down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
            down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
                       f'{gene} -og cp -out {out} -max_len {MAX_LEN} '
                       f'-stop 1 -expand 0 -rename -seq_n 10 -uniq no',
                       shell=True)
            if down.returncode == 0:
                fasta = out / 'by-gene' / f'{gene}.fasta'
                if fasta.exists:
                    yield fasta


def config(out, seed, arg):
    config = f"""Project:
-----------------------
Project name          = {out}
Type                  = chloro
Genome Range          = {arg.min}-{arg.max}
K-mer                 = {arg.kmer}
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


def clean(name):
    contigs = Path('.').glob('Contigs_*{}*'.format(name))
    options = Path('.').glob('Option_*{}*'.format(name))
    merged = Path('.').glob('Merged_contigs_{}*'.format(name))
    tmp = Path('.').glob('contigs_tmp_{}*'.format(name))
    log = Path('.').glob('log_{}.txt'.format(name))
    for i in [*contigs, *options, *merged, *tmp, *log]:
        i.rename(name/i)


def main():
    arg = parse_args()
    out = Path(Path(arg.f).stem)
    MAX_FAIL = 10
    out.mkdir()
    success = False
    fail = 0
    pattern = re.compile(r'^Assembly length\s+: +(\d+) bp$')
    for seed in get_seed(arg.taxon, out):
        log.info(f'Use {seed} as seed file.')
        config_file = config(out, seed, arg)
        test = run(f'perl NOVOPlasty2.7.2.pl -c {config_file}', shell=True)
        if test.returncode == 0:
            merged = list(Path('.').glob('Merged_contigs_{}*'.format(out)))
            if len(merged) != 0:
                with open(merged[0], 'r') as _:
                    length = re.findall(pattern, _.read())
                if len(length) != 0:
                    if min(length) >= arg.min and max(length) <= arg.max:
                        log.info(f'Assembly length (bp): {", ".join(length)}')
                        success = True
                        break
        clean(out)
        fail += 1
        if fail >= MAX_FAIL:
            break
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    with open('Failed.csv', 'a') as out:
        out.write(f'{arg.f} {arg.r}\n')
    return


if __name__ == '__main__':
    main()
