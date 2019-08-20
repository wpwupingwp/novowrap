#!/usr/bin/python3

from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from os import devnull, mkdir
from pathlib import Path
from subprocess import run
from time import sleep
from tempfile import TemporaryDirectory
from matplotlib import pyplot as plt
import argparse
import logging


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
    arg.add_argument('contig', help='contig file')
    arg.add_argument('-t', '-taxon', dest='taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
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


def down_ref(taxon, output):
    """
    Arg:
        taxon(str): given taxon name
        output(Path): output folder
    Yield:
        fasta(Path): fasta file
        out(Path): fasta file's folder
    """
    lineage = get_full_taxon(taxon)
    if lineage is not None:
        lineage = list(reversed(lineage))
    else:
        log.warning(f'Cannot find {taxon}, use Nicotiana tabacum instead.')
        lineage = list(reversed(get_full_taxon('Nicotiana tabacum')))
    if lineage is None:
        log.critical('Failed to get taxon. Quit.')
        exit(-2)
    for taxon in lineage:
        if taxon == '':
            continue
        if ' ' in taxon:
            taxon = taxon.strip('"')
            taxon = taxon.replace(' ', '_')
        # Entrez has limitation on query frenquency (3 times per second)
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
        sleep(0.5)
        query = (f'''{taxon}[Organism] AND refseq[filter] '''
                 f'''AND (chloroplast[filter] OR plastid[filter])''')
        log.info(f'Query:\t{query}')
        handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                            usehistory='y'))
        count = int(handle['Count'])
        if count == 0:
            continue
        output_file = output / f'{taxon.replace(" ", "_")}.gb'
        content = Entrez.efetch(db='nuccore', webenv=handle['WebEnv'],
                                query_key=handle['QueryKey'], rettype='gb',
                                retmode='text', retmax=1)
        output_file.write_text(content.read())
        return output_file


def blast(query, target):
    """
    Use simple BLAST with special output format.
    Args:
        query(Path): query filename
        target(Path): target filename
    Return:
        blast_out(Path): blast result filename
    """
    FMT = ('qseqid sseqid qseq sseq length pident gapopen qstart qend '
           'sstart send')
    blast_out = query.with_suffix('.blast')
    # use blastn -subject instead of makeblastdb
    blast = run(f'blastn -query {query} -subject {target} -outfmt "7 {FMT}" '
                f'-out {blast_out} -strand plus',
                shell=True, stdout=NULL, stderr=NULL)
    # remove makeblastdb result
    if blast.returncode != 0:
        log.critical('Cannot run BLAST.')
        exit(-1)
    return blast_out


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    Return [qseqid, sseqid, qseq, sseq, length, pident, gapopen, qstart, qend,
    sstart, send]
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


def get_ref_region(ref_gb, output):
    """
    Arg:
        ref_gb(Path): reference gb file, only contain one record
    Return:
        ref_region({name: [start, end]}): region location info
    """
    ref_region = {}
    for feature in SeqIO.read(ref_gb, 'gb').features:
        if (feature.type == 'misc_feature' and
                feature.qualifiers.get('software', ['', ])[0] == 'rotate_gb'):
            key = feature.qualifiers['note'][0][-4:-1]
            value = [feature.location.start, feature.location.end]
            ref_region[key] = value
#     _ = SeqIO.read(ref_gb, 'gb')
#     ref_region['All'] = [1, len(_)]

    return ref_region


def main():
    arg = parse_args()
    output = Path(Path(arg.contig).stem)
    mkdir(output)
    log.info(f'Contig:\t{arg.contig}')
    log.info(f'Taxonomy:\t{arg.taxon}')
    log.info(f'Use {output} as output folder.')

    ref_gb = down_ref(arg.taxon, output)
    run_rotate = run(f'python3 rotate_gb.py {ref_gb}', shell=True)
    if run_rotate.returncode != 0:
        exit(-2)
    new_ref_gb = list(output.glob('*.new.gb'))[0]
    ref_region_info = get_ref_region(new_ref_gb, output)
    print(ref_region_info)
    for key, value in ref_region_info.items():
        plt.plot(value, [10, 10], label=key)
    plt.show()
    ref_fasta = list(output.glob('*.rotate'))[0]
    # new gb with region info
    new_ref_gb = list(output.glob('*.new.gb'))[0]
    blast_result = blast(Path(arg.contig), ref_fasta)
    print('qseqid, sseqid, pident, gapopen, qstart, qend, sstart, send')
    for query in parse_blast_tab(blast_result):
        for i in query:
            (qseqid, sseqid, qseq, sseq, length, pident, gapopen, qstart,
             qend, sstart, send) = i
            print(qseqid, sseqid, length, pident, gapopen, qstart, qend,
                  sstart, send)
    TMP.cleanup()
    NULL.close()
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
