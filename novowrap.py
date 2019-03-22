#!/usr/bin/python3

from Bio import Entrez
from pathlib import Path
from subprocess import run
import argparse
import logging


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
    arg.add_argument('-min', default=120_000, help='minimum genome size (kb)')
    arg.add_argument('-max', default=200_000, help='minimum genome size (kb)')
    arg.add_argument('-reads_len', default=150, help='reads length')
    arg.add_argument('-seed', help='seed file')
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


def get_seq(lineage, output):
    for taxon in reversed(lineage):
        if taxon == '':
            continue
        query_str = f'{taxon}[organism] AND chloroplast[filter]'
        search = Entrez.read(Entrez.esearch(db='nuccore', term=query_str,
                                            usehistory='y'))
        if search['Count'] == '0':
            continue
        else:
            data = Entrez.efetch(db='nuccore', webenv=search['WebEnv'],
                                 query_key=search['QueryKey'], rettype='gb',
                                 retmode='text', retmax=10)
            with open(output, 'wt') as out:
                out.write(data.read())
            return True
    return False


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
    contigs = Path('.').glob('Contigs_*{}*'.format(name))
    options = Path('.').glob('Option_*{}*'.format(name))
    merged = Path('.').glob('Merged_contigs_{}*'.format(name))
    tmp = Path('.').glob('contigs_tmp_{}*'.format(name))
    log = Path('.').glob('log_{}.txt'.format(name))
    for i in [*contigs, *options, *merged, *tmp, *log]:
        i.rename(name/i)


def main():
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
    return


if __name__ == '__main__':
    main()
