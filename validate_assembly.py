#!/usr/bin/python3

from pathlib import Path
import argparse
import logging

from Bio import SeqIO
from matplotlib import pyplot as plt

from utils import down_ref, blast, parse_blast_tab, rotate_seq


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
    arg.add_argument('-n', type=int, default=5,
                     help='top n of records to keep, 0 for all')
    return arg.parse_args()


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
        alpha = 0.1
    elif old < 95:
        alpha = 0.15
    elif old < 100:
        alpha = 0.2
    else:
        alpha = 0.3
    return alpha


def draw(fasta, query, subject, ref_region, data):
    """
    Draw figure.
    Args:
        fasta(Path): fasta filename
        query(Path): query record name
        subject(Path): subject record name
        ref_region(list): reference region information
        data(list): BLAST result
    Return:
        pdf(Path): figure file
    """
    plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f'BLAST validation of {fasta}-{query} to {subject}')
    plt.xlabel('Base')
    for key, value in ref_region.items():
        plt.plot(value, [0.8, 0.8], marker='+', label=key, linewidth=10)
    plt.plot(0.5, 0.5, 'r-+', label='plus')
    plt.plot(0.5, 0.5, 'g-|', label='minus')
    plt.ylim([0.5, 1.1])
    plt.xlim(left=0)
    plt.yticks([0.65, 0.8, 0.95], labels=['minus', 'ref', 'plus'])
    plt.legend(loc='upper right')
    for i in data:
        qstart, qend, sstart, send, sstrand, pident = i
        if sstrand == 'plus':
            plt.plot([qstart, qend], [0.95, 0.95], 'r-+', linewidth=5)
            # ignore these line
            if send-sstart < 100:
                continue
            plt.fill([qstart, sstart, send, qend], [0.8, 0.95, 0.95, 0.8],
                     color='r', alpha=get_alpha(pident))
        else:
            plt.plot([qstart, qend], [0.65, 0.65], 'g-|', linewidth=5)
            plt.fill([qstart, sstart, send, qend], [0.65, 0.8, 0.8, 0.65],
                     color='#88cc88', alpha=get_alpha(pident))
    pdf = Path(f"{fasta.stem}-{Path(query+'-'+subject).with_suffix('.pdf')}")
    plt.savefig(pdf)
    plt.close()
    return pdf


def main():
    """
    Use BLAST to validate assembly result.
    Only handle first record in file.
    """
    arg = parse_args()
    arg.contig = Path(arg.contig)
    output = Path(arg.contig.stem)
    output.mkdir()
    log.info(f'Contig:\t{arg.contig}')
    log.info(f'Taxonomy:\t{arg.taxon}')
    log.info(f'Use {output} as output folder.')

    if arg.ref_gb is None:
        ref_gb = down_ref(arg.taxon)
        ref_gb.rename(output/ref_gb)
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
        r_gb, r_contig, r_regions = rotate_seq(arg.contig)
        contig_files.append(r_contig)
    if arg.n != 0:
        skip = len(contigs) - arg.n
        if skip > 0:
            log.critical(f'Skip {skip} records.')
            contig_files = contig_files[:arg.n]

    new_ref_gb, ref_fasta, ref_regions = rotate_seq(ref_gb)
    ref_region_info = get_ref_region(new_ref_gb, output)

    for i in contig_files:
        log.info(f'Analyze {i}.')
        result = compare(i, ref_fasta, output, arg)
        pdf = draw(i, result[0][0], ref_gb_name, ref_region_info, result[0][1])
        log.info(f'Write figure {pdf}.')
        # to be continued

    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
