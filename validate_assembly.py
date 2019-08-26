#!/usr/bin/python3

from Bio import Entrez, SeqIO
from os import devnull, mkdir
from pathlib import Path
from time import sleep
from tempfile import TemporaryDirectory
from matplotlib import pyplot as plt
import argparse
import logging

from utils import rotate_seq, blast, parse_blast_tab, get_full_taxon


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
    arg.add_argument('-r', '-ref_gb', dest='ref_gb', help='reference gb')
    arg.add_argument('-t', '-taxon', dest='taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
    arg.add_argument('-i', '-perc_identity', dest='perc_identity', type=float,
                     default=70.0,
                     help='minimum percentage of identity of BLAST, 0-100')
    arg.add_argument('-l', '-len_diff', dest='len_diff', type=float,
                     default=10,
                     help='maximum percentage of length differnce of query to'
                     'reference, 0-100')
    return arg.parse_args()


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


def compare(query, reference, output, arg):
    results = []
    blast_result = blast(Path(query), reference, arg.perc_identity)
    # only one record in file, loop only for unpack
    for query in parse_blast_tab(blast_result):
        record = []
        for i in query:
            (qseqid, sseqid, qseq, sseq, sstrand, length, pident, gapopen,
             qstart, qend, sstart, send) = i
            record.append([qstart, qend, sstart, send, sstrand, pident])
            # print(qseqid, sseqid, length, pident, gapopen, qstart, qend,
            #       sstart, send)
        results.append([qseqid, record])
    return results


def get_alpha(old):
    """
    Given 0-100, return 0, 0.5, 0.75, 0.95, 1
    """
    alpha = 0
    if old < 50:
        alpha = 0
    elif old < 80:
        alpha = 0.05
    elif old < 95:
        alpha = 0.1
    elif old < 100:
        alpha = 0.15
    else:
        alpha = 0.2
    return alpha


def draw(contig, query, subject, ref_region, data):
    """
    Draw figure.
    """
    plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f'BLAST validation of {contig}-{query} to {subject}')
    plt.xlabel('Base')
    for key, value in ref_region.items():
        plt.plot(value, [0.8, 0.8], marker='+', label=key, linewidth=10)
    plt.plot(0.5, 0.5, 'r-+', label='plus')
    plt.plot(0.5, 0.5, 'g-|', label='minus')
    plt.ylim([0.5, 1.1])
    plt.xlim(left=0)
    plt.yticks([0.7, 0.8, 0.9], label=['minus', 'ref', 'plus'])
    plt.legend(loc='upper right')
    for i in data:
        qstart, qend, sstart, send, sstrand, pident = i
        if sstrand == 'plus':
            plt.plot([qstart, qend], [0.9, 0.9], 'r-+', linewidth=5)
            plt.fill_between([min(qstart, sstart), max(qend, send)],
                             [0.8, 0.8], [0.9, 0.9],
                             alpha=get_alpha(pident),
                             color='#ff8888')
        else:
            plt.plot([qstart, qend], [0.7, 0.7], 'g-|', linewidth=5)
            plt.fill_between([min(qstart, sstart), max(qend, send)],
                             [0.7, 0.7], [0.8, 0.8],
                             alpha=get_alpha(pident),
                             color='#88cc88')
    plt.savefig(f"{contig.stem}-{Path(query+'-'+subject).with_suffix('.pdf')}")
    plt.close()


def main():
    """
    Use BLAST to validate assembly result.
    Only handle first record in file.
    """
    arg = parse_args()
    arg.contig = Path(arg.contig)
    output = Path(arg.contig.stem)
    mkdir(output)
    log.info(f'Contig:\t{arg.contig}')
    log.info(f'Taxonomy:\t{arg.taxon}')
    log.info(f'Use {output} as output folder.')

    if arg.ref_gb is None:
        ref_gb = down_ref(arg.taxon, output)
    else:
        ref_gb = arg.ref_gb
    _ = SeqIO.read(ref_gb, 'gb')
    ref_gb_name = _.name
    # make folder clean
    with open(output / (ref_gb_name+'.gb'), 'w') as d, open(ref_gb, 'r') as s:
        d.write(s.read())
    ref_gb = output / (ref_gb_name + '.gb')
    ref_len = len(SeqIO.read(ref_gb, 'gb'))

    contigs = list(SeqIO.parse(arg.contig, 'fasta'))
    contig_files = []
    if len(contigs) > 1:
        log.warning(f'Find {len(contigs)} records in {arg.contig}.')
        log.info('Divide them into different files.')
        for idx, record in enumerate(contigs):
            filename = arg.contig.with_suffix(f'.{idx}.fasta')
            record_len = len(record)
            if abs(1-(record_len/ref_len))*100 > arg.len_diff:
                log.warning(f'The length difference of record with reference'
                            f'({abs(record_len-ref_len)} bp) is out of limit'
                            f'({arg.len_diff}%).')
                new_filename = str(filename) + '.bad_length'
                SeqIO.write(record, new_filename, 'fasta')
                log.warning(f'Skip {new_filename}.')
                continue
            log.info(f'\t{filename}')
            SeqIO.write(record, filename, 'fasta')
            r_gb, r_contig, r_regions = rotate_seq(filename)
            r_gb.rename(output/r_gb)
            r_regions.rename(output/r_regions)
            contig_files.append(r_contig)
    else:
        contig_files.append(arg.contig)

    new_ref_gb, ref_fasta, ref_regions = rotate_seq(ref_gb)
    ref_region_info = get_ref_region(new_ref_gb, output)

    for i in contig_files:
        result = compare(i, ref_fasta, output, arg)
        draw(i, result[0][0], ref_gb_name, ref_region_info, result[0][1])

    TMP.cleanup()
    NULL.close()
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
