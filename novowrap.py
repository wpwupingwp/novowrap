#!/usr/bin/python3

from Bio import Entrez, SeqIO
from os import devnull, remove
from pathlib import Path
from platform import system
from subprocess import run
from tempfile import TemporaryDirectory
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
# temporary directory
TMP = TemporaryDirectory()


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


def get_seq(taxon, output, gene=None):
    """
    Use BarcodeFinder to get seed or reference sequence.
    Arg:
        taxon(str): given taxon name
        output(Path): output folder
        gene(tuple): gene name
    Return:
        fasta(Path): fasta file
        out(Path): fasta file's folder
    """
    if gene is None:
        genes = ('rbcL', 'matK', 'psaB', 'psaC', 'rrn23')
    else:
        genes = (gene, )
    # +, -, -, -, +
    MAX_LEN = 250000
    if system() == 'Windows':
        python = 'python'
    else:
        python = 'python3'
    last_taxon = ''
    null_handle = open(devnull, 'w')
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
                           shell=True, stdout=null_handle, stderr=null_handle)
            else:
                down = run(f'{python} -m BarcodeFinder -taxon {taxon} -gene '
                           f'{gene} -og cp -out {out} -max_len {MAX_LEN} '
                           f'-stop 1 -expand 0 -rename -seq_n 10 -uniq no ',
                           shell=True, stdout=null_handle, stderr=null_handle)
            if down.returncode == 0:
                fasta = out / 'by-gene' / f'{gene}.fasta'
                if fasta.exists:
                    yield fasta, out
            else:
                log.critical('Failed to run BarcodeFinder.')
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


def repeat_and_reverse(fasta, taxon):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.
    Detect direction of sequence by BLASTing rbcL, reverse if needed.
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
                log.warning(f'Detected reversed sequence {sseqid}.')
    new_fasta = fasta.with_suffix('.new')
    new = []
    for i in SeqIO.parse(fasta, 'fasta'):
        i.seq = i.seq + i.seq
        # negative strand or not found by blast
        if strand.get(i.id, '-') == '-':
            i.seq = i.seq[::-1]
            log.warning(f'Posible negative strand found: {i.id}.')
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
    # check blast
    with open(devnull, 'w') as out:
        _ = run(f'makeblastdb -in {target} -dbtype nucl -out {target}',
                shell=True, stdout=out)
    if _.returncode != 0:
        exit('Cannot run BLAST.')
    blast_out = query.with_suffix('.blast')
    blast = run(f'blastn -query {query} -db {target} -outfmt "7 {FMT}" -out '
                f'{blast_out}', shell=True)
    # remove makeblastdb result
    remove(str(target)+'.nhr')
    remove(str(target)+'.nin')
    remove(str(target)+'.nsq')
    if blast.returncode != 0:
        log.critical('Cannot run BLAST.')
    return blast_out


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    Return [qseqid, sseqid, qseq, sseq, pident, gapopen, qstart, qend, sstart,
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
        yield query


def rotate(fasta, taxon):
    """
    Rotate sequences, from LSC (trnH-psbA) to IRa, SSC, IRb.
    Arg:
        fasta(Path or str): fasta filename
        taxon(str): fasta's taxon
    Return:
        new_fasta???

    """
    # FMT = 'qseqid sseqid qseq sseq pident gapopen qstart qend sstart send'
    if not isinstance(fasta, Path):
        fasta = Path(fasta)
    repeat_fasta = repeat_and_reverse(fasta, taxon)
    blast_result = blast(repeat_fasta, repeat_fasta)
    for query in parse_blast_tab(blast_result):
        max_aln_len = 0
        locations = []
        name = ''
        length = 0
        seqs = []
        for hit in query:
            # print(hit[-4:])
            (qseqid, sseqid, qseq, sseq, qlen, pident, gapopen,
             qstart, qend, sstart, send) = hit
            name = qseqid
            seqs.append(qseq)
            length = raw_qlen = qlen // 2
            # mismatch
            if pident != 100 or qseqid != sseqid or gapopen != 0:
                continue
            location = tuple(sorted([qstart, qend, sstart, send]))
            # hit across origin and repeat
            if location[-1] - location[0] > raw_qlen:
                continue
            aln_len = abs(qstart-qend) + 1
            # short hit or hit of whole sequence or whole repeat sequence
            if aln_len < max_aln_len or aln_len in (qlen, raw_qlen):
                continue
            else:
                max_aln_len = aln_len
            if location not in locations:
                locations.append(location)
            else:
                continue
        if not locations:
            continue
        locations.sort(key=lambda x: x[0])
        # remove extra hit, first two if enough
        locations = locations[:2]
        print('\n{}: {} {}\t{} {}\t{}'.format(name, *locations[0], length))
        print('{}: {} {}\t{} {}\t{}'.format(name, *locations[1], length))
        # ira_start, ira_end, irb_start, irb_end
        a = slice(locations[0][0]-1, locations[0][1])
        b = slice(locations[0][1], locations[0][2]-1)
        c = slice(locations[0][2]-1, locations[0][3])
        d = slice(locations[1][1], locations[1][2]-1)
        if (locations[0][2]-locations[0][1]) > (
                locations[1][2]-locations[1][1]):
            region_LSC = b
            region_IRa = c
            region_SSC = d
            region_IRb = a
        else:
            region_LSC = d
            region_IRa = a
            region_SSC = b
            region_IRb = c
        print(f'{name}, LSC {region_LSC}, IRa {region_IRa}, SSC {region_SSC},'
              f'IRb {region_IRb}')
        old_seq = sorted(seqs, key=lambda x: len(x))[-1]
        new_seq = ''.join([old_seq[region_LSC], old_seq[region_IRa],
                           old_seq[region_SSC], old_seq[region_IRb]])
        # remove repeat
        half_old_seq = old_seq[0:len(old_seq)//2]
        print('old', len(half_old_seq), 'new', len(new_seq))
        with open('/tmp/test.fasta', 'w') as out:
            out.write('>old\n{}\n'.format(half_old_seq))
            out.write('>new\n{}\n'.format(new_seq))
            out.write(f'>lsc\n{old_seq[region_LSC]}\n')
            out.write(f'>ssc\n{old_seq[region_SSC]}\n')
            out.write(f'>ira\n{old_seq[region_IRa]}\n')
            out.write(f'>irb\n{old_seq[region_IRb]}\n')

    remove(repeat_fasta)
    remove(blast_result)
    # return new


def neaten_out(source, dest):
    """
    Organize NOVOPlasty output.
    Return sequence lengths and fasta filename list.
    According to seed, adjust direction of assembled sequences.
    Arg:
        source(Path): current directory
        dest(Path): directory to move
    Return:
        seq_len(list): sequences length
        fasta(list): assembled fasta
    """
    def merge_to_fasta(merge):
        options = []
        fasta = merge.with_suffix('.fasta')
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

    contigs = list(source.glob('Contigs_*'))
    options = list(source.glob('Option_*'))
    merged = list(source.glob('Merged_contigs_*'))
    print(merged)
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
    return seq_len, fasta


def main():
    arg = parse_args()
    out = Path(Path(arg.f).stem)
    out.mkdir()
    success = False
    fail = 0
    rbcL_list = []
    for seed, folder in get_seq(arg.taxon, out):
        # collect rbcL
        if str(seed).endswith('rbcL.fasta'):
            rbcL_list.append(seed)
        log.info(f'Use {seed} as seed file.')
        config_file = config(out, seed, arg)
        run(f'perl NOVOPlasty2.7.2.pl -c {config_file}', shell=True)
        # novoplasty generates outputs in current folder
        # use rbcL to detect strand direction
        seq_len, assembled = neaten_out(Path().cwd(), folder)
        if len(seq_len) != 0:
            # f-string cannot use *
            log.info('Assembled length:\t{}.'.format(*seq_len))
            if min(seq_len) >= arg.min and max(seq_len) <= arg.max:
                rotated, rbcL_list = rotate(assembled, arg.taxon, rbcL_list)
                # break
    # rotate
    #             break
        else:
            log.warn('Assembled failed.')
        fail += 1
        if fail >= arg.try_n:
            break
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    with open('Failed.csv', 'a') as out:
        out.write(f'{arg.f} {arg.r}\n')
    log.info(f'Clean temporary files {TMP.name}.')
    TMP.cleanup()
    return


if __name__ == '__main__':
    main()
