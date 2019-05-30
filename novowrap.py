#!/usr/bin/python3

from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
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
TEMP_LOG = 'novowrap.log'
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
NULL = open(devnull, 'w')


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
    # strand: +, -, -, -, +
    candidate_genes = ('rbcL', 'matK', 'psaB', 'psaC', 'rrn23')
    if gene is None:
        genes = candidate_genes
    else:
        genes = [gene, ]
        genes.extend(candidate_genes)
    log.info(f'Candidate seed genes: {", ".join(genes)}.')
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
            log.info(f'Querying {gene} of {taxon}.')
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
            else:
                log.critical('Failed to run BarcodeFinder. Retry...')
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
    # # check blast
    # _ = run(f'makeblastdb -in {target} -dbtype nucl -out {target}',
    #         shell=True, stdout=NULL, stderr=NULL)
    # if _.returncode != 0:
    #     # if BLAST fail, directly quit
    #     log.critical('Cannot run BLAST.')
    #     log.info('Exit.')
    #     exit(-1)
    blast_out = query.with_suffix('.blast')
    blast = run(f'blastn -query {query} -subject {target} -outfmt "7 {FMT}" '
                f'-out {blast_out}', shell=True)
    # remove makeblastdb result
    #remove(str(target)+'.nhr')
    #remove(str(target)+'.nin')
    #remove(str(target)+'.nsq')
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
        success(bool): success or not
        new_fasta(Path): rotated file name
    """
    # FMT = 'qseqid sseqid qseq sseq pident gapopen qstart qend sstart send'
    if not isinstance(fasta, Path):
        fasta = Path(fasta)
    repeat_fasta = repeat_and_reverse(fasta, taxon)
    blast_result = blast(repeat_fasta, repeat_fasta)
    # analyze blast result
    new_fasta = fasta.with_suffix('.rotate')
    new_regions = fasta.with_suffix('.regions')
    new_gb = fasta.with_suffix('.gb')
    for query in parse_blast_tab(blast_result):
        max_aln_len = 0
        locations = []
        name = ''
        seqs = []
        for hit in query:
            (qseqid, sseqid, qseq, sseq, qlen, pident, gapopen,
             qstart, qend, sstart, send) = hit
            name = qseqid
            seqs.append(qseq)
            raw_qlen = qlen // 2
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
        if len(locations) < 2:
            remove(repeat_fasta)
            remove(blast_result)
            return fasta, False
        # remove extra hit, first two if enough
        locations = locations[:2]
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
        old_seq = sorted(seqs, key=lambda x: len(x))[-1]
        seq_LSC = old_seq[region_LSC]
        seq_IRa = old_seq[region_IRa]
        seq_SSC = old_seq[region_SSC]
        seq_IRb = old_seq[region_IRb]
        log.info(f'{name}, LSC {region_LSC}, IRa {region_IRa}, '
                 f'SSC {region_SSC}, IRb {region_IRb}')
        new_seq = ''.join([seq_LSC, seq_IRa, seq_SSC, seq_IRb])
        # remove repeat
        half_old_seq = old_seq[0:len(old_seq)//2]
        log.info(f'old, {len(half_old_seq)}, new, {len(new_seq)}')
        record = SeqRecord(Seq(new_seq, alphabet=IUPAC.ambiguous_dna))
        record.annotations['accession'] = 'Unknown'
        record.annotations['organism'] = name
        # output
        offset = -1
        for f_name, f in zip(('LSC', 'IRa', 'SSC', 'IRb'),
                             (region_LSC, region_IRa, region_SSC, region_IRb)):
            length = f.stop - f.start
            record.features.append(SeqFeature(
                FeatureLocation(offset+1, length+offset+1),
                type='misc_feature',
                qualifiers={'name': f_name, 'note': f_name}, strand=1))
            offset += length
        assert seq_SSC == str(record.features[2].extract(record).seq)
        with open(new_fasta, 'a') as out:
            # out.write('>old\n{}\n'.formak(half_old_seq))
            out.write(f'>{name}\n{new_seq}\n')
        with open(new_regions, 'a') as out2:
            out2.write(f'>{name}-LSC\n{seq_LSC}\n')
            out2.write(f'>{name}-IRa\n{seq_IRa}\n')
            out2.write(f'>{name}-SSC\n{seq_SSC}\n')
            out2.write(f'>{name}-IRb\n{seq_IRb}\n')
        with open(new_gb, 'a') as out3:
            SeqIO.write(record, out3, 'gb')
        with open('validate.f', 'w') as out:
            for i in record.features:
                SeqIO.write(i.extract(record), out, 'fasta')

    remove(repeat_fasta)
    remove(blast_result)
    return new_fasta, True


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
        fasta = merge.with_suffix('.long.fasta')
        with open(merge, 'r') as raw:
            for line in raw:
                if line.startswith('>'):
                    options.append([line, next(raw)])
        if options:
            with open(fasta, 'w') as out:
                for i in options:
                    out.write(i[0])
                    out.write(i[1])
            return Path(fasta), max([len(i[1]) for i in options])
        else:
            return None, 0

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
    # circularized file use different file format with Merged file -_-
    circularized = list(source.glob('Circularized_assembly*'))
    for i in circularized:
        # assume each circularized file only have one fasta record
        with open(i, 'r') as _:
            raw = _.readlines()
            raw = [i for i in raw if not i.startswith('>')]
            n = sum([len(i) for i in raw])
        seq_len.append(n)
    fasta.extend(circularized)
    tmp = list(source.glob('contigs_tmp_*'))
    log = list(source.glob('log_*.txt'))
    for i in (*contigs, *options, *merged, *tmp, *log):
        i.replace(dest/i.name)
    assembled = []
    for i in fasta:
        new_loc = dest / i.name
        i.replace(new_loc)
        assembled.append(new_loc)
    return seq_len, assembled


def main():
    log.info('Welcome to novowrap.')
    arg = parse_args()
    out = Path(Path(arg.f).stem+'-out').absolute()
    log.info(f'Use {out} as output folder.')
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
        log.info(f'Organize NOVOPlasty output of {seed}.')
        seq_len, assembled = neaten_out(Path().cwd(), folder)
        if len(seq_len) != 0:
            # f-string cannot use *
            log.info('Assembled length:\t{}.'.format(*seq_len))
            if min(seq_len) < arg.min or max(seq_len) > arg.max:
                log.warning('Length does not fulfill requirement.')
            else:
                rotate_result = [rotate(i, arg.taxon) for i in assembled]
                # (filename, True/False)
                rotated = [i for i in rotate_result if i[1] is True]
                if rotated:
                    success = True
                    break
                else:
                    log.warning('Cannot find correct conformation.')
        else:
            log.warn(f'Assembled with {seed} failed.')
        fail += 1
        if fail >= arg.try_n:
            break
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    with open('Failed.csv', 'a') as out:
        out.write(f'{arg.f} {arg.r}\n')
    log.info(f'Clean temporary files.')
    TMP.cleanup()
    NULL.close()
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
