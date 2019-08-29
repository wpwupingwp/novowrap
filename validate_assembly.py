#!/usr/bin/python3

from pathlib import Path
import argparse
import logging

from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np

from utils import down_ref, blast, parse_blast_tab, rotate_seq, rc_region


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
                     default=0.7,
                     help='minimum percentage of identity of BLAST, 0-100')
    arg.add_argument('-l', '-len_diff', dest='len_diff', type=float,
                     default=0.1,
                     help='maximum percentage of length differnce of query to'
                     'reference, 0-100')
    arg.add_argument('-n', dest='top', type=int, default=0,
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
        perc_identity(float): percent identity for BLAST, need multiply 100
    Return:
        results[0][1]: BLAST data
    """
    results = []
    blast_result = blast(Path(query), reference, perc_identity*100)
    # only one record in file, loop is for unpack
    for query in parse_blast_tab(blast_result):
        record = []
        for i in query:
            (qseqid, sseqid, qseq, sseq, sstrand, length, pident, gapopen,
             qstart, qend, sstart, send) = i
            record.append([qstart, qend, sstart, send, sstrand, pident])
        results.append([qseqid, record])
    assert len(results) == 1
    return results[0][1]


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


def draw(title, ref_region, option_region, data):
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
    plt.title(f"The validation result of {title.replace('-', ' and ')}")
    plt.xlabel('Base')
    for key, value in ref_region.items():
        plt.plot(value[:2], [0.8, 0.8], marker='+', label=key, linewidth=10)
    for key, value in option_region.items():
        plt.plot([value[1], value[1]], [0.93, 0.97],
                 'k--', linewidth=2, alpha=0.3)
        plt.plot([value[1], value[1]], [0.63, 0.67],
                 'k--', linewidth=2, alpha=0.3)
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
            plt.fill([sstart, qstart, qend, send], [0.8, 0.95, 0.95, 0.8],
                     color='r', alpha=get_alpha(pident))
        # minus strand
        else:
            plt.plot([qstart, qend], [0.65, 0.65], 'g-|', linewidth=5)
            plt.fill([send, sstart, qend, qstart], [0.8, 0.8, 0.65, 0.65],
                     color='#88cc88', alpha=get_alpha(pident))
    pdf = Path(title).with_suffix('.pdf')
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


def divide_records(fasta, output, ref_len, len_diff=0.1, top=0):
    """
    Make sure each file has only one record.
    Args:
        fasta(Path or str): fasta file
        output(Path): output folder
        ref_len(int): length of reference, to filter bad records
        len_diff: maximum allowed length difference
        top(int): number of records to keep
    Returns:
        option_files(list(Path)):  list of divided files
    """
    options = list(SeqIO.parse(fasta, 'fasta'))
    fasta = Path(fasta)
    option_files = []
    if len(options) > 1:
        log.warning(f'Found {len(options)} records in {fasta}.')
        log.info('Divide them into different files.')
        for idx, record in enumerate(options):
            filename = output / f'{idx}_{fasta}'
            record_len = len(record)
            if abs(1-(record_len/ref_len)) > len_diff:
                log.critical(f'The length difference of record with reference '
                             f'({abs(record_len-ref_len)} bp) is out of limit '
                             f'({len_diff:.0%}, {len_diff*ref_len:.0f}).')
                new_filename = str(filename) + '.bad_length'
                SeqIO.write(record, new_filename, 'fasta')
                log.warning(f'Skip {new_filename}.')
                continue
            SeqIO.write(record, filename, 'fasta')
            option_files.append(rotate_seq(filename))
    else:
        option_files.append(clean_rotate(fasta, output))
    if top != 0:
        skip = len(option_files) - top
        if skip > 0:
            log.critical(f'Skip {skip} records.')
            option_files = option_files[:top]
    return option_files


def main():
    """
    Use BLAST to validate assembly result.
    Only handle first record in file.
    """
    arg = parse_args()
    arg.contig = Path(arg.contig)
    output = Path(arg.contig.stem)
    output.mkdir()
    validated = []
    result_info = []
    log.info(f'Contig:\t{arg.contig}')
    log.info(f'Taxonomy:\t{arg.taxon}')
    log.info(f'Use {output} as output folder.')

    if arg.ref_gb is None:
        ref_gb = down_ref(arg.taxon, output)
    else:
        ref_gb = output / Path(arg.ref_gb)
        # Path.rename is problematic
        with open(arg.ref_gb, 'r') as i, open(ref_gb, 'w') as o:
            o.write(i.read())
    ref_len = len(SeqIO.read(ref_gb, 'gb'))

    option_files = divide_records(arg.contig, output, ref_len, arg.len_diff,
                                  arg.top)

    # ref already in output
    new_ref_gb, ref_fasta, ref_regions = rotate_seq(ref_gb)
    ref_region_info = get_region(new_ref_gb)

    for i in option_files:
        i_gb, i_fasta, i_regions = i
        log.info(f'Analyze {i_fasta}.')
        option_region_info = get_region(i_gb)
        option_len = len(SeqIO.read(i_fasta, 'fasta'))
        compare_result = compare(i_fasta, ref_fasta, arg.perc_identity)
        fig_title = str(output / f'{i_fasta.stem}-{ref_gb.stem}')
        pdf = draw(fig_title, ref_region_info, option_region_info,
                   compare_result)
        log.info(f'Write figure {pdf}.')
        # to be continued
        log.info('Detecting reverse complement region.')

        plus = np.zeros(option_len, dtype=bool)
        minus = np.zeros(option_len, dtype=bool)
        count = {}
        for region in option_region_info:
            count[region] = {'loc': slice(*option_region_info[region][:2]),
                             'len': option_region_info[region][2],
                             'plus': 0, 'minus': 0}

        for hsp in compare_result:
            qstart, qend, sstart, send, sstrand, pident = hsp
            if sstrand == 'plus':
                plus[qstart-1:qend] = True
            else:
                minus[qstart-1:qend] = True
        # count bases
        for rgn in count:
            p_slice = plus[count[rgn]['loc']]
            m_slice = minus[count[rgn]['loc']]
            count[rgn]['plus'] = np.count_nonzero(p_slice)
            count[rgn]['minus'] = np.count_nonzero(m_slice)
            count[rgn]['union'] = np.count_nonzero(p_slice | m_slice)
        # assign strand
        for rgn in count:
            # also use arg.perc_identity for region
            min_rgn_len = count[rgn]['len'] * arg.perc_identity
            p = count[rgn]['plus']
            m = count[rgn]['minus']
            u = count[rgn]['union']
            if p == m == 0:
                count[rgn]['strand'] = 'missing'
                log.critical(f'Region {rgn} of {i_fasta} is missing.')
            elif u < min_rgn_len:
                count[rgn]['strand'] = 'incomplete'
                log.critical(f'Region {rgn} of {i_fasta} is incomplete.')
            elif p > min_rgn_len:
                count[rgn]['strand'] = 'plus'
            elif m > min_rgn_len:
                count[rgn]['strand'] = 'minus'
        # do not rc IR
        if count['LSC']['strand'] == count['SSC']['strand'] == 'minus':
            log.info(f'Reverse complement the whole sequence of {i_fasta}.')
            edited = rc_region(i_fasta, i_regions, 'whole')
            pass
        elif count['LSC']['strand'] == 'minus':
            log.warning(f'Reverse complement the LSC of {i_fasta}.')
            edited = rc_region(i_fasta, i_regions, 'LSC')
            pass
        elif count['SSC']['strand'] == 'minus':
            log.warning(f'Reverse complement the SSC of {i_fasta}.')
            edited = rc_region(i_fasta, i_regions, 'SSC')
        else:
            edited = i_fasta
        if edited != i_fasta:
            new_compare_result = compare(edited, ref_fasta, arg.perc_identity)
            fig_title = str(output / f'{edited.stem}-{ref_gb.stem}')
            pdf = draw(fig_title, ref_region_info, option_region_info,
                       new_compare_result)
            log.info(f'Write figure {pdf}.')
        validated.append(edited)

    log.info('Validated sequences:')
    for i in validated:
        log.info(f'\t{i}')
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
