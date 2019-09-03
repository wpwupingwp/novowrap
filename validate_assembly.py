#!/usr/bin/python3

from pathlib import Path
import argparse
import logging

from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np

from utils import down_ref, blast, parse_blast_tab
from utils import get_fmt, rotate_seq, get_regions, rc_regions

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
    arg.add_argument('input', help='input filename')
    arg.add_argument('-r', '-ref', dest='ref', help='reference gb')
    arg.add_argument('-t', '-taxon', dest='taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
    arg.add_argument('-i', '-perc_identity', dest='perc_identity', type=float,
                     default=0.7,
                     help='minimum percentage of identity of BLAST, 0-100')
    arg.add_argument('-l', '-len_diff', dest='len_diff', type=float,
                     default=0.2,
                     help='maximum percentage of length differnce of query to'
                     'reference, 0-100')
    return arg.parse_args()


def move(source, dest):
    """
    Move source to dest and return dest.
    Args:
        source(Path): old path
        dest(Path or str): new path
    Return:
        dest(Path): new path
    """
    source.rename(dest)
    return Path(dest)


def divide_records(fasta, output, ref_len, len_diff=0.1):
    """
    Make sure each file has only one record.
    Args:
        fasta(Path or str): fasta file
        output(Path): output folder
        ref_len(int): length of reference, to filter bad records
        len_diff: maximum allowed length difference
    Returns:
        option_files(list(Path)):  list of divided files
        info(list): file info
    """
    options = list(SeqIO.parse(fasta, 'fasta'))
    fasta = Path(fasta)
    divided = {}
    keys = ('skip,gb,fasta,length,LSC,IRa,SSC,IRb,missing,incomplete,rc,'
            'figure,figure_after').split(',')
    log.warning(f'Found {len(options)} records in {fasta}.')
    log.info('Divide them into different files.')
    for idx, record in enumerate(options):
        skip = False
        r_gb = r_fasta = ''
        filename = output / f'{idx}_{fasta}'
        divided[filename] = dict((key, '') for key in keys)
        record_len = len(record)
        record_len_diff = abs(1-(record_len/ref_len))
        if record_len_diff > len_diff:
            log.critical(f'The length difference of NO.{idx+1} record '
                         f'with reference is out of limit '
                         f'({record_len_diff:.2%} > {len_diff:.2%}).')
            skip = True
        SeqIO.write(record, filename, 'fasta')
        if not skip:
            r_gb, r_fasta = rotate_seq(filename)
            if r_gb is not None:
                divided[filename].update({'gb': r_gb, 'fasta': r_fasta,
                                          'length': record_len, 'skip': skip})
            else:
                divided[filename]['skip'] = True
    return divided


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
            (qseqid, sseqid, sstrand, length, pident, gapopen, qstart, qend,
             sstart, send) = i
            record.append([qstart, qend, sstart, send, sstrand, pident])
        results.append(record)
    assert len(results) == 1
    blast_result.unlink()
    return results[0]


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


def draw(title, ref_regions, option_regions, data):
    """
    Draw figure.
    Args:
        title(str): figure title
        ref_regions(dict): reference region information
        option_regions(dict): option sequence's region information
        data(list): BLAST result
    Return:
        pdf(Path): figure file
    """
    ignore_offset = len(ref_regions['IRa'])*2 + len(ref_regions['SSC'])
    plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f"The validation result of {title.replace('-', ' and ')}",
              pad=10)
    plt.xlabel('Base')
    for key, value in ref_regions.items():
        plt.plot([value.location.start, value.location.end], [0.8, 0.8],
                 marker='+', label=key, linewidth=10)
    for key, value in option_regions.items():
        plt.plot([value.location.end, value.location.end], [0.93, 0.97],
                 'k--', linewidth=2, alpha=0.3)
        plt.plot([value.location.end, value.location.end], [0.63, 0.67],
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
            # ignore these line
            if abs(qstart-sstart) > ignore_offset:
                continue
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


def validate_regions(length, regions, compare, perc_identity=0.7):
    """
    Use BLAST results to validate regions.
    Args:
        length(int): length of whole sequences
        regions(dict): regions information
        compare(list): BLAST result
        perc_identity: threshold of region similarity
    Returns:
        count(dict): count information
        to_rc(str or None): regions need to reverse-complement
    """
    count = {}
    to_rc = None
    plus = np.zeros(length, dtype=bool)
    minus = np.zeros(length, dtype=bool)
    for region in regions:
        count[region] = {'plus': 0, 'minus': 0}
    for hsp in compare:
        qstart, qend, sstart, send, sstrand, pident = hsp
        if sstrand == 'plus':
            plus[qstart-1:qend] = True
        else:
            minus[qstart-1:qend] = True
    # count bases
    for rgn in count:
        start = int(regions[rgn].location.start)
        end = int(regions[rgn].location.end)
        p_slice = plus[start:end]
        m_slice = minus[start:end]
        count[rgn]['plus'] = np.count_nonzero(p_slice)
        count[rgn]['minus'] = np.count_nonzero(m_slice)
        count[rgn]['union'] = np.count_nonzero(p_slice | m_slice)
    # assign strand
    for rgn in count:
        # also use arg.perc_identity for region
        min_rgn_len = len(regions[rgn]) * perc_identity
        p = count[rgn]['plus']
        m = count[rgn]['minus']
        u = count[rgn]['union']
        if p == m == 0:
            count[rgn]['strand'] = 'missing'
        elif u < min_rgn_len:
            count[rgn]['strand'] = 'incomplete'
        elif p > min_rgn_len:
            count[rgn]['strand'] = 'plus'
        elif m > min_rgn_len:
            count[rgn]['strand'] = 'minus'
    # do not rc IR
    to_rc = None
    if count['LSC']['strand'] == count['SSC']['strand'] == 'minus':
        to_rc = 'whole'
    elif count['LSC']['strand'] == 'minus':
        to_rc = 'LSC'
    elif count['SSC']['strand'] == 'minus':
        to_rc = 'SSC'
    return count, to_rc


def main():
    """
    Use BLAST to validate assembly result.
    """
    arg = parse_args()
    arg.input = Path(arg.input)
    output = Path(arg.input.stem)
    output.mkdir()
    tmp = output / 'Temporary'
    tmp.mkdir()
    log.info(f'Contig:\t{arg.input}')
    if arg.ref is not None:
        log.info(f'Reference:\t{arg.ref}')
        arg.taxon = None
    else:
        log.info(f'Taxonomy:\t{arg.taxon}')
    log.info(f'Use {output} as output folder.')
    # get ref
    if arg.ref is None:
        arg.taxon, ref_gb = down_ref(arg.taxon, output)
        fmt = 'gb'
    else:
        fmt = get_fmt(arg.ref)
        ref_gb = move(Path(arg.ref), output/arg.ref)
    ref_len = len(SeqIO.read(ref_gb, fmt))
    # ref already in output
    new_ref_gb, ref_fasta = rotate_seq(ref_gb)
    if new_ref_gb is None:
        log.critical('Cannot get rotated reference sequence.')
        log.critical('Please consider to use another reference.')
        raise SystemExit
    ref_regions = get_regions(new_ref_gb)
    divided = divide_records(arg.input, output, ref_len, arg.len_diff)
    for i in divided:
        success = False
        divided[i]['success'] = success
        if divided[i]['skip']:
            continue
        i_gb = divided[i]['gb']
        i_fasta = divided[i]['fasta']
        log.info(f'Analyze {i_fasta}.')
        option_regions = get_regions(i_gb)
        # add regions info
        for _ in option_regions:
            divided[i][_] = len(option_regions[_])
        compare_result = compare(i_fasta, ref_fasta, arg.perc_identity)
        fig_title = str(output / f'{i_fasta.stem}-{ref_gb.stem}')
        pdf = draw(fig_title, ref_regions, option_regions,
                   compare_result)
        divided[i]['figure'] = pdf
        log.info('Detecting reverse complement region.')
        option_len = divided[i]['length']
        count, to_rc = validate_regions(option_len, option_regions,
                                        compare_result, arg.perc_identity)
        for rgn in count:
            if count[rgn]['strand'] == 'missing':
                divided[i]['missing'] = ','.join([divided[i]['missing'], rgn])
                log.critical(f'Region {rgn} of {i_fasta} is missing.')
                success = False
            elif count[rgn]['strand'] == 'incomplete':
                log.critical(f'Region {rgn} of {i_fasta} is incomplete.')
                divided[i]['incomplete'] = ','.join([divided[i]['incomplete'],
                                                     rgn])
                success = False
            else:
                pass
        if to_rc is not None:
            log.warning(f'Reverse complement the {to_rc} of {i_fasta}.')
            rc_gb, rc_fasta = rc_regions(i_gb, to_rc)
            # clean old files
            i_fasta = move(i_fasta, tmp/(i_fasta.name+'.tmp'))
            i_gb = move(i_gb, tmp/(i_gb.name+'.tmp'))
            rc_gb = move(rc_gb, rc_gb.parent/rc_gb.name.replace('_RC_', ''))
            rc_fasta = move(rc_fasta, rc_fasta.parent/rc_fasta.name.replace(
                '_RC_', ''))
            new_compare_result = compare(rc_fasta, ref_fasta,
                                         arg.perc_identity)
            fig_title = str(output / f'_RC_{rc_fasta.stem}-{ref_gb.stem}')
            new_regions = get_regions(rc_gb)
            pdf = draw(fig_title, ref_regions, new_regions,
                       new_compare_result)
            divided[i]['figure_after'] = pdf
            divided[i]['rc'] = to_rc
            divided[i]['gb'] = rc_gb
            divided[i]['fasta'] = rc_fasta
            for _ in new_regions:
                divided[i][_] = len(new_regions[_])
            # validate again
            count_2, to_rc_2 = validate_regions(option_len, new_regions,
                                                new_compare_result,
                                                arg.perc_identity)
            if to_rc_2 is None:
                success = True
        else:
            divided[i]['figure_after'] = divided[i]['figure']
        divided[i]['success'] = success

    log.info('Validated sequences:')
    for i in divided:
        if divided[i]['success']:
            log.info(f"\t{divided[i]['fasta']}")
    reference_info = output / f'{output.name}-Reference.csv'
    output_info = output / f'{output.name}-Results.csv'
    with open(output_info, 'w') as out:
        out.write('Raw,Success,Skip,gb,fasta,Length,LSC,IRa,'
                  'SSC,IRb,Missing,Incomplete,RC_region,Figure,Figure_after\n')
        for record in divided:
            # format is easier than f-string for dict
            simple = divided[record]
            if simple['gb'] != '':
                simple['gb'] = simple['gb'].name
                simple['fasta'] = simple['fasta'].name
                simple['figure'] = simple['figure'].name
                simple['figure_after'] = simple['figure_after'].name
            record = record.name
            out.write('{},'.format(record))
            out.write('{success},{skip},{gb},{fasta},{length},{LSC},{IRa},'
                      '{SSC},{IRb},{missing},{incomplete},{rc},{figure},'
                      '{figure_after}\n'.format(**simple))
    with open(reference_info, 'w') as out:
        out.write('Reference,Taxon,Length,LSC,IRa,SSC,IRb\n')
        out.write('{},{},{},{},{},{},{}\n'.format(
            ref_fasta.name, arg.taxon, ref_len, len(ref_regions['LSC']),
            len(ref_regions['IRa']), len(ref_regions['SSC']),
            len(ref_regions['IRb'])))
    log.info(f'Validation result was written into {output_info}')
    log.info(f'Reference information was written into {reference_info}')
    log.info('Bye.')
    return


if __name__ == '__main__':
    main()
