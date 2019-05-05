#!/usr/bin/python3

from Bio import Entrez, SeqIO
from os import devnull
from pathlib import Path
from platform import system
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
    arg.add_argument('-min', default=40000, help='minimum genome size (kb)')
    arg.add_argument('-max', default=300000, help='maximum genome size (kb)')
    arg.add_argument('-reads_len', default=150, help='reads length')
    arg.add_argument('-taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
    arg.add_argument('-try', dest='try_n', type=int, default=12,
                     help='maximum tried times')
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
    # species name contains genus
    if lineage[-1] != '' and lineage[-2] != '':
        lineage[-1] = f'"{lineage[-2]} {lineage[-1]}"'
    return lineage


def get_seed(taxon, output):
    genes = ('rbcL', 'matK', 'psaB', 'psaC', 'rrn23')
    MAX_LEN = 250000
    if system() == 'Windows':
        python = 'python'
    else:
        python = 'python3'
    last_taxon = ''
    for taxon in reversed(get_full_taxon(taxon)):
        if taxon == '':
            continue
        for gene in genes:
            if ' ' in taxon:
                taxon = taxon.strip('"')
                taxon = taxon.replace(' ', '_')
            out = output / f'{taxon}-{gene}'
            log.info(f'Querying {out}.')
            if last_taxon != '':
                down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
                           f'{gene} -og cp -out {out} -max_len {MAX_LEN} '
                           f'-stop 1 -expand 0 -rename -seq_n 10 -uniq no '
                           f'-exclude "{last_taxon}"[organism]',
                           shell=True)
            else:
                down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
                           f'{gene} -og cp -out {out} -max_len {MAX_LEN} '
                           f'-stop 1 -expand 0 -rename -seq_n 10 -uniq no ',
                           shell=True)
            if down.returncode == 0:
                fasta = out / 'by-gene' / f'{gene}.fasta'
                if fasta.exists:
                    yield fasta, out
        last_taxon = taxon


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


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    """
    query = []
    with open(filename, 'r', encoding='utf-8') as raw:
        for line in raw:
            if line.startswith('# BLAST'):
                yield query
                query.clear()
            elif line.startswith('#'):
                pass
            else:
                query.append(line.strip())


def repeat(fasta):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.
    """
    new_fasta = fasta + '.new'
    new = []
    for i in SeqIO.parse(fasta, 'fasta'):
        i.seq = i.seq + i.seq
        new.apend(i)
    SeqIO.write(new, new_fasta, 'fasta')
    return new_fasta


def rotate(fasta):
    MIN_IR = 1000
    MAX_IR = 100000
    with open(devnull, 'w') as out:
        _ = run(f'makeblastdb -in {fasta} -dbtype nucl', shell=True,
                stdout=out)
    if _.returncode != 0:
        exit('Cannot run BLAST.')
    fmt = 'qseqid sseqid qseq sseq pident gapopen qstart qend sstart send'
    blast_out = 'blast.txt'
    blast = run(f'blastn -query {fasta} -db {fasta} -outfmt "7 {fmt}" -out '
                f'{blast_out}', shell=True)
    if blast.returncode != 0:
        exit('Cannot run BLAST.')
    for query in parse_blast_tab(blast_out):
        for hit in query:
            (qseqid, sseqid, qseq, sseq, pident, gapopen,
             qstart, qend, sstart, send, *_) = hit.split('\t')
            if float(pident) != 100 or qseqid != sseqid or int(gapopen) != 0:
                continue
            qlen = abs(int(qstart)-int(qend)) + 1
            if qlen > MAX_IR or qlen < MIN_IR:
                continue
            else:
                print(qseqid, sseqid, qlen, qstart, qend, sstart, send)
                print('next')
                # break

    # return new


def merge_to_fasta(merge):
    options = []
    fasta = str(merge) + '.fasta'
    with open(merge, 'r') as raw:
        for line in raw:
            if line.startswith('>'):
                options.append([line, next(raw)])
    if options:
        with open(fasta, 'w') as out:
            for i in options:
                out.write(i[0])
                out.write(i[1])
        return Path(fasta), max([len(i) for i in options])
    else:
        return None, 0


def clean(source, dest):
    contigs = list(source.glob('Contigs_*'))
    options = list(source.glob('Option_*'))
    merged = list(source.glob('Merged_contigs_*'))
    fasta = []
    seq_len = []
    for i in merged:
        f, s = merge_to_fasta(i)
        if f is not None:
            fasta.append(f)
            seq_len.append(s)
    tmp = list(source.glob('contigs_tmp_*'))
    log = list(source.glob('log_*.txt'))
    for i in [*contigs, *options, *merged, *tmp, *log, *fasta]:
        i.rename(dest/i)
    return seq_len


def main():
    arg = parse_args()
    out = Path(Path(arg.f).stem)
    out.mkdir()
    success = False
    fail = 0
    for seed, folder in get_seed(arg.taxon, out):
        log.info(f'Use {seed} as seed file.')
        config_file = config(out, seed, arg)
        run(f'perl NOVOPlasty2.7.2.pl -c {config_file}', shell=True)
        seq_len = clean(Path('.'), folder)
        if len(seq_len) != 0:
            # f-string cannot use *
            log.info('Assembled length:\t{}.'.format(*seq_len))
            if min(seq_len) >= arg.min and max(seq_len) <= arg.max:
                break
        else:
            log.warn('Assembled failed.')
        fail += 1
        if fail >= arg.try_n:
            break
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    with open('Failed.csv', 'a') as out:
        out.write(f'{arg.f} {arg.r}\n')
    return


if __name__ == '__main__':
    main()
