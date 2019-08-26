#!/usr/bin/python3

from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from os import devnull
from pathlib import Path
from platform import system
from subprocess import run
from time import sleep
from tempfile import TemporaryDirectory
import argparse
import logging

from rotate import rotate_seq


# temporary directory
TMP = TemporaryDirectory()
NULL = open(devnull, 'w')
# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger(__name__)
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('-f', required=True, help='forward fastq file')
    arg.add_argument('-r', required=True, help='reverse fastq file')
    arg.add_argument('-kmer', choices=range(23, 40, 2), default=39, type=int,
                     help='kmer size')
    arg.add_argument('-min', default=100000, help='minimum genome size (KB)')
    arg.add_argument('-max', default=200000, help='maximum genome size (KB)')
    arg.add_argument('-mem', default=30, type=int, help='maximum memory (GB)')
    arg.add_argument('-reads_len', default=150, help='reads length')
    arg.add_argument('-taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
    arg.add_argument('-try', dest='try_n', type=int, default=5,
                     help='maximum tried times')
    # arg.add_argument('-split', default=1_000_000,
    #                  help='reads to use (million), set to 0 to skip split')
    return arg.parse_args()


def get_full_taxon(taxon):
    """
    Get full lineage of given taxon
    Return lineage list, only contains Kingdom, Phylum, Class, Order, Family,
    Genus, Species.
    Arg:
        taxon(str): given taxon name, could be common name, use quotation mark
        if necessary
    Return:
        lineage(list): lineage list, from higher to lower rank
    """
    split = taxon.split(' ')
    if len(split) >= 2 and split[0][0].isupper() and split[0][1].islower():
        name = ' '.join(split[0:2])
    else:
        name = split[0]
    Entrez.email = 'guest@example.org'
    search = Entrez.read(Entrez.esearch(db='taxonomy', term=f'"{name}"'))
    if search['Count'] == '0':
        if ' ' not in name:
            log.critical(f'Cannot find {name} in NCBI Taxonomy.')
            return None
        if ' ' in name:
            name = split[0]
            sleep(0.5)
            search = Entrez.read(Entrez.esearch(db='taxonomy',
                                                term=f'"{name}"'))
            if search['Count'] == '0':
                log.critical(f'Cannot find {name} in NCBI Taxonomy.')
                return None
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


def get_seq(taxon, output, gene=None):
    """
    Use BarcodeFinder to get seed or reference sequence.
    Arg:
        taxon(str): given taxon name
        output(Path): output folder
        gene(tuple): gene name
    Yield:
        fasta(Path): fasta file
        out(Path): fasta file's folder
    """
    # strand: +, -, -, -, +
    candidate_genes = ('rbcL', 'matK', 'psaB', 'psaC', 'rrn23')
    lineage = get_full_taxon(taxon)
    if lineage is not None:
        lineage = list(reversed(lineage))
    else:
        log.warning(f'Cannot find {taxon}, use Nicotiana tabacum instead.')
        lineage = list(reversed(get_full_taxon('Nicotiana tabacum')))
    if lineage is None:
        log.critical('Failed to get taxon. Quit.')
        exit(-2)
    if gene is None:
        genes = candidate_genes
    else:
        genes = [gene, ]
        genes.extend(candidate_genes)
    if system() == 'Windows':
        python = 'python'
    else:
        python = 'python3'
    for gene in genes:
        last_taxon = ''
        for taxon in lineage:
            if taxon == '':
                continue
            if ' ' in taxon:
                taxon = taxon.strip('"')
                taxon = taxon.replace(' ', '_')
            out = output / f'{gene}-{taxon}'
            # Entrez has limitation on query frenquency (3 times per second)
            # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
            sleep(0.5)
            if last_taxon != '':
                down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
                           f'{gene} -og cp -out {out} -stop 1 -expand 0 '
                           f'-rename -seq_n 10 -uniq no '
                           f'-exclude "{last_taxon}"[organism]',
                           shell=True, stdout=NULL, stderr=NULL)
            else:
                down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
                           f'{gene} -og cp -out {out} '
                           f'-stop 1 -expand 0 -rename -seq_n 10 -uniq no ',
                           shell=True, stdout=NULL, stderr=NULL)
            if down.returncode == 0:
                fasta = out / 'by-gene' / f'{gene}.fasta'
                if fasta.exists:
                    yield fasta, out
                    # if nearest seed fail, higher rank may be useless, too
                    break
            else:
                log.warning(f'Cannot get {gene} of {taxon}. Retry...')
            last_taxon = taxon


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


def repeat_and_reverse(fasta, taxon):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.
    Detect direction of sequence by BLASTing rbcL, reverse if needed.
    Replace illegal character in sequence with "N".
    Args:
        fasta(Path): fasta filename
        taxon(Path): taxonomy of given fasta
    Return:
        new_fasta(Path): new_fasta's name
    """
    strand = {}
    # get rbcL reference
    rbcL = Path(TMP.name) / f'{taxon}-rbcL_ref' / 'by-gene' / 'rbcL.fasta'
    if rbcL.exists():
        pass
    else:
        rbcL, _ = next(get_seq(taxon, Path(TMP.name), 'rbcL'))

    rbcL_blast = blast(rbcL, fasta)
    for query in parse_blast_tab(rbcL_blast):
        for hit in query:
            (qseqid, sseqid, qseq, sseq, qlen, pident, gapopen,
             qstart, qend, sstart, send) = hit
            if qstart < qend and sstart < send:
                strand[sseqid] = '+'
            else:
                strand[sseqid] = '-'
    new_fasta = fasta.with_suffix('.new')
    new = []
    for i in SeqIO.parse(fasta, 'fasta'):
        n_star = str(i.seq).count('*')
        if n_star != 0:
            log.info(f'Remove {n_star} "*" in {i.description}.')
            i.seq = Seq(str(i.seq).replace('*', ''))
        i.seq = i.seq + i.seq
        # negative strand or not found by blast
        if strand.get(i.id, '-') == '-':
            # reverse complement instead of reverse
            old_name = i.description
            i = i.reverse_complement()
            i.id = old_name
            i.description = ''
            log.warning(f'Found reversed sequence in "{i.id}". Reverse back.')
        new.append(i)
    SeqIO.write(new, new_fasta, 'fasta')
    return new_fasta


def blast(query, target):
    """
    Use simple BLAST with special output format.
    Args:
        query(Path): query filename
        target(Path): target filename
    Return:
        blast_out(Path): blast result filename
    """
    FMT = 'qseqid sseqid qseq sseq qlen pident gapopen qstart qend sstart send'
    blast_out = query.with_suffix('.blast')
    # use blastn -subject instead of makeblastdb
    blast = run(f'blastn -query {query} -subject {target} -outfmt "7 {FMT}" '
                f'-out {blast_out}', shell=True, stdout=NULL, stderr=NULL)
    # remove makeblastdb result
    if blast.returncode != 0:
        log.critical('Cannot run BLAST.')
        exit(-1)
    return blast_out


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    Return [qseqid, sseqid, qseq, sseq, qlen, pident, gapopen, qstart, qend, sstart,
    send]
    Arg:
        filename(Path): blast result file
    Return:
        line(list): parsed result
    """
    query = []
    with open(filename, 'r', encoding='utf-8') as raw:
        for line in raw:
            if line.startswith('# BLAST'):
                if query:
                    yield query
                query.clear()
            elif line.startswith('#'):
                pass
            else:
                line = line.strip().split('\t')
                line[4:] = list(map(float, line[4:]))
                line[4:] = list(map(int, line[4:]))
                query.append(line)


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
            out.write(line)
    return new


def neaten_out(source, dest):
    """
    Organize NOVOPlasty output.
    Return sequence lengths and fasta filename list.
    According to seed, adjust direction of assembled sequences.
    Arg:
        source(Path): current directory
        dest(Path): directory to move
    Return:
        assembled(list): assembled fasta
    """
    contigs = list(source.glob('Contigs_*'))
    options = list(source.glob('Option_*'))
    merged = list(source.glob('Merged_contigs_*'))
    circularized = list(source.glob('Circularized_assembly*'))
    tmp = list(source.glob('contigs_tmp_*'))
    log = list(source.glob('log_*.txt'))
    assembled = []
    # move to dest folder
    for i in (*options, *contigs, *tmp, *log):
        i.replace(dest/i.name)
    # move to dest folder, generate clean fasta
    for i in (*merged, *circularized):
        i.replace(dest/i.name)
        new_loc = dest / i.name
        fasta = txt_to_fasta(new_loc)
        assembled.append(fasta)
    return assembled


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

    success = False
    fail = 0
    for seed, folder in get_seq(arg.taxon, out):
        if fail >= arg.try_n:
            log.critical(f'Too much failure. Quit.')
            break
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
        assembled = neaten_out(Path().cwd(), folder)
        if len(assembled) == 0:
            log.warning(f'Assembled with {seed.name} failed.')
            fail += 1
            continue
        rotate_result = [rotate_seq(i) for i in assembled]
        if any(rotate_result):
            success = True
            break
        else:
            log.warning('Cannot find correct conformation for all assembly of '
                        f'{seed.name}.')
            fail += 1
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    TMP.cleanup()
    NULL.close()
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
