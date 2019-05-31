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
            log.warning(f'Replace {n_star} illegal character in '
                        f'{i.description} with "N".')
            i.seq = Seq(str(i.seq).replace('*', 'N'))
        i.seq = i.seq + i.seq
        # negative strand or not found by blast
        if strand.get(i.id, '-') == '-':
            i.seq = i.seq[::-1]
            log.warning(f'Detected reversed sequence {i.name}. Reverse back.')
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


def rotate(fasta, taxon, min_len=40000, max_len=300000):
    """
    Rotate sequences, from LSC (trnH-psbA) to IRa, SSC, IRb.
    Repeat length parameters for easily call rotate function instead of run
    whole program.
    Arg:
        fasta(Path or str): fasta filename
        taxon(str): taxon of fasta
        min_len: minimum length of sequence
        max_len: maximum length of sequence
    Return:
        success(bool): success or not
    """
    # FMT = 'qseqid sseqid qseq sseq pident gapopen qstart qend sstart send'
    if not isinstance(fasta, Path):
        fasta = Path(fasta)
    seq_len = [len(i) for i in SeqIO.parse(fasta, 'fasta')]
    if not seq_len or min(seq_len) < min_len or max(seq_len) > max_len:
        log.warning("The sequences' length of assembly {fasta}"
                    "failed to meet requirement. Skip.")
        return False
    repeat_fasta = repeat_and_reverse(fasta, taxon)
    blast_result = blast(repeat_fasta, repeat_fasta)
    # analyze blast result
    new_fasta = fasta.with_suffix('.rotate')
    new_regions = fasta.with_suffix('.regions')
    new_gb = fasta.with_suffix('.gb')
    success = False
    # blast and fasta use same order
    for query, seq in zip(parse_blast_tab(blast_result),
                          SeqIO.parse(repeat_fasta, 'fasta')):
        locations = set()
        # use name to get uniform sequence id with BLAST
        # biopython 's seq.name and seq.description is complex
        name = ''
        original_seq_len = len(seq) // 2
        ambiguous_base_n = len(str(seq).strip('ATCGatcg')) // 2
        # only use hit of IR to IR
        max_aln_n = 0
        # because of ambiguous base, percent of identity bases may not be 100%
        if ambiguous_base_n != 0:
            log.warning(f'Found {ambiguous_base_n} ambiguous bases.')
            p_ident_min = int((1-(ambiguous_base_n/original_seq_len))*100)
        else:
            p_ident_min = 100
        for hit in query:
            (qseqid, sseqid, qseq, sseq, qlen, pident, gapopen,
             qstart, qend, sstart, send) = hit
            name = qseqid
            # only self
            if qseqid != sseqid or gapopen != 0:
                continue
            # mismatch
            location = tuple(sorted([qstart, qend, sstart, send]))
            if pident < p_ident_min or qseqid != sseqid or gapopen != 0:
                continue
            # hit across origin and repeat
            if location[-1] - location[0] > original_seq_len:
                continue
            aln_len = abs(qstart-qend) + 1
            # self to self or self to repeat self
            if len(set(location)) != 4:
                continue
            # origin to repeat
            if aln_len == len(seq):
                continue
            # filter short hit
            if aln_len < max_aln_n:
                continue
            else:
                max_aln_n = aln_len
            locations.add(location)
        if not locations:
            continue
        locations = list(locations)
        locations.sort(key=lambda x: x[1]-x[0], reverse=True)
        locations.sort(key=lambda x: x[0])
        if len(locations) < 2:
            remove(repeat_fasta)
            remove(blast_result)
            return False
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
        # self to self is the longest match
        seq_LSC = seq[region_LSC]
        seq_IRa = seq[region_IRa]
        seq_SSC = seq[region_SSC]
        seq_IRb = seq[region_IRb]
        log.info(f'Original region of {name}:')
        log.info(f'LSC {region_LSC}, IRa {region_IRa}, '
                 f'SSC {region_SSC}, IRb {region_IRb}')
        new_seq = seq_LSC + seq_IRa + seq_SSC + seq_IRb
        new_seq.seq.alphabet = IUPAC.ambiguous_dna
        assert len(seq)//2 == len(new_seq)
        new_seq.annotations['accession'] = 'Unknown'
        new_seq.annotations['organism'] = name
        # output
        offset = -1
        for f_name, f in zip(('LSC', 'IRa', 'SSC', 'IRb'),
                             (region_LSC, region_IRa, region_SSC, region_IRb)):
            length = f.stop - f.start
            new_seq.features.append(SeqFeature(
                FeatureLocation(offset+1, length+offset+1),
                type='misc_feature',
                qualifiers={'name': f_name}, strand=1))
            offset += length
        assert str(seq_SSC.seq) == str(
            new_seq.features[2].extract(new_seq).seq)
        with open(new_fasta, 'a') as out:
            out.write(f'>{name}\n{new_seq}\n')
        with open(new_regions, 'a') as out2:
            out2.write(f'>{name}-LSC\n{seq_LSC}\n')
            out2.write(f'>{name}-IRa\n{seq_IRa}\n')
            out2.write(f'>{name}-SSC\n{seq_SSC}\n')
            out2.write(f'>{name}-IRb\n{seq_IRb}\n')
        with open(new_gb, 'a') as out3:
            SeqIO.write(new_seq, out3, 'gb')
        success = True

    remove(repeat_fasta)
    remove(blast_result)
    if not success:
        return False
    log.info(f'Rotated {fasta.name} to uniform conformation {new_fasta.name}.')
    log.info(f'Chloroplast region information were written into {new_gb}')
    return True


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
    for i in (*contigs, *tmp, *log):
        i.replace(dest/i.name)
    # move to dest folder, generate clean fasta
    for i in (*options, *merged, *circularized):
        new_loc = dest / i.name
        new_name = new_loc.with_suffix('.convert')
        i.replace(new_loc)
        SeqIO.convert(new_loc, 'fasta', new_name, 'fasta')
        assembled.append(new_name)
    return assembled


def main():
    log.info('Welcome to novowrap.')
    arg = parse_args()
    out = Path(Path(arg.f).stem+'-out').absolute()
    log.info(f'Use {out} as output folder.')
    out.mkdir()
    success = False
    fail = 0
    for seed, folder in get_seq(arg.taxon, out):
        log.info(f'Use {seed} as seed file.')
        config_file = config(out, seed, arg)
        run_novo = run(f'perl NOVOPlasty2.7.2.pl -c {config_file}', shell=True)
        if run_novo.returncode != 0:
            log.critical('Failed to run NOVOPlasty. Quit.')
            exit(-1)
        # novoplasty generates outputs in current folder
        # use rbcL to detect strand direction
        log.info(f'Organize NOVOPlasty output of {seed.name}.')
        # novoplasty use current folder as output folder
        assembled = neaten_out(Path().cwd(), folder)
        if len(assembled) == 0:
            log.warn(f'Assembled with {seed} failed.')
            continue
        rotate_result = [rotate(i, arg.taxon) for i in assembled]
        if any(rotate_result):
            success = True
            break
        else:
            log.warning('Cannot find correct conformation for all assembly of '
                        f'{seed}.')
        fail += 1
        if fail >= arg.try_n:
            log.critical(f'Too much failure, quit.')
            break
    if not success:
        log.critical(f'Failed to assemble {arg.f} and {arg.r}.')
    TMP.cleanup()
    NULL.close()
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
