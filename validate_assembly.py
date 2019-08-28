#!/usr/bin/python3

from pathlib import Path
import argparse
import logging

from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np

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


def get_region(gb):
    """
    Arg:
        gb(Path): rotate_seq generated gb file, only contains one record
    Return:
        region({name: [start, end, length]}): region location info
    """
    ref_region = {}
    for feature in SeqIO.read(gb, 'gb').features:
        if (feature.type == 'misc_feature' and
                feature.qualifiers.get('software', ['', ])[0] == 'rotate_seq'):
            key = feature.qualifiers['note'][0][-4:-1]
            value = [feature.location.start, feature.location.end,
                     len(feature)]
            ref_region[key] = value
    return ref_region


def compare(query, reference, perc_identity):
    """
    Use BLAST to compare two records.
    Args:
        query(Path or str): query file
        reference(Path or str): reference file
        arg(Path): BLAST parameters
    Return:
        results[0][0]: qseqid
        results[0][1]: BLAST data
    """
    results = []
    blast_result = blast(Path(query), reference, perc_identity)
    # only one record in file, loop is for unpack
    for query in parse_blast_tab(blast_result):
        record = []
        for i in query:
            (qseqid, sseqid, qseq, sseq, sstrand, length, pident, gapopen,
             qstart, qend, sstart, send) = i
            record.append([qstart, qend, sstart, send, sstrand, pident])
        results.append([qseqid, record])
    assert len(results) == 1
    return results[0][0], results[0][1]


def get_alpha(old):
    """
    Given 0-100, return 0, 0.5, 0.75, 0.95, 1
    Args:
        old(float): percent of identity
    Return:
        alpha(float): alpha value
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


def draw(title, ref_region, data):
    """
    Draw figure.
    Args:
        title(str): figure title
        ref_region(list): reference region information
        data(list): BLAST result
    Return:
        pdf(Path): figure file
    """
    ignore_offset = (ref_region['IRa'][1] - ref_region['IRa'][0])*2 + (
        ref_region['SSC'][1]-ref_region['SSC'][0])
    plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f'BLAST validation of {title}')
    plt.xlabel('Base')
    for key, value in ref_region.items():
        plt.plot(value[:2], [0.8, 0.8], marker='+', label=key, linewidth=10)
    # no repeat legend
    plt.plot(0.5, 0.5, 'r-+', linewidth=5, label='Plus')
    plt.plot(0.5, 0.5, 'g-|', linewidth=5, label='Minus')
    plt.ylim([0.5, 1.1])
    plt.xlim(left=0)
    plt.yticks([0.65, 0.8, 0.95], labels=['Minus', 'Reference', 'Plus'])
    plt.legend(loc='upper right')
    for i in data:
        qstart, qend, sstart, send, sstrand, pident = i
        if sstrand == 'plus':
            plt.plot([qstart, qend], [0.95, 0.95], 'r-+', linewidth=5)
            if abs(qstart-sstart) > ignore_offset:
                continue
            # ignore these line
            plt.fill([qstart, sstart, send, qend], [0.8, 0.95, 0.95, 0.8],
                     color='r', alpha=get_alpha(pident))
        else:
            plt.plot([qstart, qend], [0.65, 0.65], 'g-|', linewidth=5)
            plt.fill([qstart, send, sstart, qend], [0.65, 0.8, 0.8, 0.65],
                     color='#88cc88', alpha=get_alpha(pident))
    title = Path(title)
    pdf = title.with_suffix('.pdf')
    plt.savefig(pdf)
    plt.close()
    return pdf


def clean_rotate(filename, output):
    """
    Make the folder clean.
    """
    r_gb, r_contig, r_regions = rotate_seq(filename)
    for i in r_gb, r_contig, r_regions:
        # rename does not return new name
        i.rename(output/i)
    return output/r_gb, output/r_contig, output/r_regions


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
        dest = output / ref_gb
        ref_gb.rename(dest)
        ref_gb = dest
    else:
        ref_gb = arg.ref_gb
    _ = SeqIO.read(ref_gb, 'gb')
    ref_gb_name = _.name
    # make folder clean
    with open(output / (ref_gb_name+'.gb'), 'w') as d, open(ref_gb, 'r') as s:
        d.write(s.read())
    ref_gb = output / (ref_gb_name + '.gb')
    ref_len = len(SeqIO.read(ref_gb, 'gb'))

    options = list(SeqIO.parse(arg.contig, 'fasta'))
    option_files = []
    if len(options) > 1:
        log.warning(f'Find {len(options)} records in {arg.contig}.')
        log.info('Divide them into different files.')
        for idx, record in enumerate(options):
            filename = output / f'{idx}-{arg.contig}'
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
            option_files.append(rotate_seq(filename))
    else:
        option_files.append(clean_rotate(arg.contig, output))
    if arg.n != 0:
        skip = len(option_files) - arg.n
        if skip > 0:
            log.critical(f'Skip {skip} records.')
            option_files = option_files[:arg.n]

    # ref already in output
    new_ref_gb, ref_fasta, ref_regions = rotate_seq(ref_gb)
    ref_region_info = get_region(new_ref_gb)

    for i in option_files:
        i_gb, i_fasta, i_regions = i
        log.info(f'Analyze {i_fasta}.')
        qseqid, compare_result = compare(i_fasta, ref_fasta, arg.perc_identity)
        fig_title = output / f'{i_fasta.stem}_{qseqid}-{ref_gb_name}'
        pdf = draw(fig_title, ref_region_info, compare_result)
        log.info(f'Write figure {pdf}.')
        # to be continued
        plus = {}
        minus = {}
        log.info('Detecting reverse complement region.')
        for region in ref_region_info:
            plus[region] = np.zeros(ref_region_info[region][2], dtype=bool)
            minus[region] = np.zeros(ref_region_info[region][2], dtype=bool)
        for hsp in compare_result:
            qstart, qend, sstart, send, sstrand, pident = hsp

        # np.count_nonzero(lsc_plus[20:30])

    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
