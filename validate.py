#!/usr/bin/python3

from os import environ
from pathlib import Path
import argparse
import logging

from Bio import SeqIO
from matplotlib import use as mpl_use
if environ.get('DISPLAY', '') == '':
    mpl_use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from utils import get_ref, blast, parse_blast_tab, move
from utils import get_fmt, rotate_seq, get_regions, rc_regions


# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
# inherit logger from novowrap, if not called by it, doesn't matter to
# name the logger 'novowrap'
log = logging.getLogger('novowrap')


def parse_args(arg_list=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('input', help='input filename')
    arg.add_argument('-r', '-ref', dest='ref', help='reference gb')
    arg.add_argument('-t', '-taxon', dest='taxon', default='Nicotiana tabacum',
                     help='Taxonomy name')
    options = arg.add_argument_group('Option')
    options.add_argument('-i', '-perc_identity', dest='perc_identity',
                         type=float, default=0.7,
                         help='minimum percentage of identity of BLAST, 0-100')
    options.add_argument('-l', '-len_diff', dest='len_diff', type=float,
                         default=0.2, help='maximum percentage of length '
                         'differnce of query to' 'reference, 0-100')
    options.add_argument('-s', '-seed', dest='seed',
                         help='seed used in assembly, only for caller')
    options.add_argument('-o', '-out', dest='out', help='output folder')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def divide_records(fasta, output, ref_len, tmp, len_diff=0.1):
    """
    Make sure each file has only one record.
    Args:
        fasta(Path or str): fasta file
        output(Path): output folder
        ref_len(int): length of reference, to filter bad records
        tmp(Path): temp folder
        len_diff: maximum allowed length difference
    Returns:
        option_files(list(Path)):  list of divided files
        info(list): file info
    """
    def _insert_suffix(old, suffix):
        """
        Insert suffix before old's suffix
        old.old_suffix -> old.suffix.old_suffix
        """
        old_suffix = old.suffix
        new = old.with_suffix(suffix)
        new = Path(str(new)+old_suffix)
        return new

    options = list(SeqIO.parse(fasta, 'fasta'))
    divided = {}
    keys = ('gb,fasta,length,LSC,IRa,SSC,IRb,missing,incomplete,'
            'rc,figure,figure_after,skip').split(',')
    log.info(f'Found {len(options)} records in {fasta.name}.')
    if len(options) > 1:
        log.info('Divide them into different files.')
    log.info(f"Check record's length (reference: {ref_len} bp, "
             f"difference limit {len_diff:.2%}).")
    for idx, record in enumerate(options):
        skip = False
        if len(options) > 1:
            filename = output / _insert_suffix(fasta, f'.{idx+1}').name
        else:
            filename = output / fasta.name
        divided[filename] = dict((key, '') for key in keys)
        record_len = len(record)
        record_len_diff = (record_len/ref_len) - 1
        divided[filename]['fasta'] = filename
        divided[filename]['length'] = record_len
        divided[filename]['length_diff'] = record_len_diff
        if abs(record_len_diff) > len_diff:
            log.warning(f'Skip NO.{idx+1} record ({record_len} bp, '
                        f'length difference {record_len_diff:.2%}).')
            skip = 'undersize' if record_len_diff < 0 else 'oversize'
        SeqIO.write(record, filename, 'fasta')
        if not skip:
            r_gb, r_fasta = rotate_seq(filename, tmp=tmp)
            if r_gb is not None:
                divided[filename].update({'gb': r_gb, 'fasta': r_fasta,
                                          'length': record_len})
                move(filename,
                     tmp/filename.with_suffix('.raw').name)
            else:
                skip = 'structure_unusual'
        divided[filename]['skip'] = skip
        log.debug(f'{skip}')
    return divided


def compare(query, reference, tmp, perc_identity):
    """
    Use BLAST to compare two records.
    Args:
        query(Path or str): query file
        reference(Path or str): reference file
        tmp(Path): temp folder
        perc_identity(float): percent identity for BLAST, need multiply 100
    Return:
        results[0][1]: BLAST data
    """
    results = []
    blast_result, blast_log = blast(Path(query), reference, perc_identity*100)
    if blast_result is None:
        return None
    # only one record in file, loop is for unpack
    for query in parse_blast_tab(blast_result):
        record = []
        for i in query:
            (qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen,
             qstart, qend, sstart, send) = i
            record.append([qstart, qend, sstart, send, sstrand, pident])
        results.append(record)
    # assert len(results) == 1
    move(blast_result, tmp/blast_result)
    move(blast_log, tmp/blast_log)
    return results[0]


def get_alpha(old):
    """
    Given 0-100, return 0, 0.5, 0.75, 0.95, 1
    Args:
        old(float): percent of identity
    Return:
        alpha(float): alpha value
    """
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


def draw(ref_gb, seq_gb, data):
    """
    Draw figure.
    Args:
        ref_gb(Path): reference genbank file
        seq_gb(Path): sequence genbank file
        data(list): BLAST result
    Return:
        pdf(Path): figure file
    """
    ref_regions = get_regions(ref_gb)
    seq_regions = get_regions(seq_gb)
    title = f'{seq_gb.stem} and {ref_gb.stem}'
    ignore_offset = len(ref_regions['IRa'])*2 + len(ref_regions['SSC'])
    plt.rcParams.update({'font.size': 16, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f'Validation of {title}', pad=10)
    plt.xlabel('Base')
    for key, value in ref_regions.items():
        plt.plot([value.location.start, value.location.end], [0.8, 0.8],
                 marker='+', label=key, linewidth=10)
    for key, value in seq_regions.items():
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
    # pdf = seq_gb.with_name(seq_gb.stem+'.pdf')
    pdf = seq_gb.with_suffix('.pdf')
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
        strand_info(dict): region's strand information
    """
    count = {}
    bad_region = False
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
        elif p >= min_rgn_len:
            count[rgn]['strand'] = 'plus'
        elif m >= min_rgn_len:
            count[rgn]['strand'] = 'minus'
        else:
            count[rgn]['strand'] = 'unknown'
    # do not rc IR
    to_rc = None
    if count['LSC']['strand'] == count['SSC']['strand'] == 'minus':
        to_rc = 'whole'
    elif count['LSC']['strand'] == 'minus':
        to_rc = 'LSC'
    elif count['SSC']['strand'] == 'minus':
        to_rc = 'SSC'
    strand_info = {}
    for rgn in count:
        if count[rgn]['strand'] == 'missing':
            if strand_info.get('missing', '') == '':
                strand_info['missing'] = rgn
            else:
                strand_info['missing'] += f'-{rgn}'
            log.critical(f'Region {rgn} is missing.')
            bad_region = True
        elif count[rgn]['strand'] == 'incomplete':
            if strand_info.get('incomplete', '') == '':
                strand_info['incomplete'] = rgn
            else:
                strand_info['incomplete'] += f'-{rgn}'
            log.critical(f'Region {rgn} is incomplete.')
            bad_region = True
        else:
            pass
    return count, to_rc, strand_info, bad_region


def validate_main(arg_str=None):
    """
    Use BLAST to validate assembly result.
    Args:
        arg_str(str): arguments string
    Return:
        validated(list): list contains validated rotated fasta files
        output_info(str): result csv, empty string for failed result
    """
    validated = []
    output_info = ''

    if arg_str is None:
        arg = parse_args()
    else:
        arg = parse_args(arg_str.split(' '))
    arg.input = Path(arg.input).absolute()
    if arg.out is None:
        output = Path(arg.input.stem+'-out').absolute()
    else:
        output = Path(arg.out).absolute()
    if not output.exists():
        output.mkdir()
    tmp = output / 'Temp'
    if not tmp.exists():
        tmp.mkdir()
    log.info(f'Input:\t{arg.input}')
    if arg.ref is not None:
        log.info(f'Reference:\t{arg.ref}')
        fmt = get_fmt(arg.ref)
        ref_gb = Path(arg.ref)
        ref_gb = move(ref_gb, tmp/ref_gb.name, copy=True)
    else:
        log.info(f'Taxonomy:\t{arg.taxon}')
        ref_gb, ref_taxon = get_ref(arg.taxon, tmp)
        if ref_gb is None:
            log.critical('Failed to get reference.')
            log.debug(f'{arg.input} {arg.ref} REF_NOT_FOUND\n')
            return validated, output_info
        ref_gb = move(ref_gb, tmp/ref_gb.name)
        fmt = 'gb'
    log.debug(f'Use {output} as output folder.')
    ref_len = len(SeqIO.read(ref_gb, fmt))
    r_ref_gb, r_ref_fasta = rotate_seq(ref_gb, tmp=tmp)
    if r_ref_gb is None:
        return validated, output_info
    ref_regions = get_regions(r_ref_gb)
    if r_ref_gb is None:
        log.critical('Cannot get rotated reference sequence.')
        log.critical('Please consider to use another reference.')
        log.debug(f'{arg.input} {arg.ref} REF_CANNOT_ROTATE\n')
        return validated, output_info
    divided = divide_records(arg.input, output, ref_len, tmp, arg.len_diff)
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
        compare_result = compare(i_fasta, r_ref_fasta, tmp, arg.perc_identity)
        if compare_result is None:
            log.critical('Cannot run BLAST.')
            log.debug(f'{arg.input} {arg.ref} BLAST_FAIL\n')
            return validated, output_info
        pdf = draw(r_ref_gb, i_gb, compare_result)
        pdf = move(pdf, output/pdf.name)
        log.info('Detecting reverse complement region.')
        option_len = divided[i]['length']
        count, to_rc, strand_info, bad_region = validate_regions(
            option_len, option_regions, compare_result, arg.perc_identity)
        divided[i].update(strand_info)
        if bad_region:
            # skip sequences with bad region
            continue
        if to_rc is not None:
            log.warning(f'Reverse complement the {to_rc} of {i_fasta.name}.')
            rc_fasta = rc_regions(i_gb, to_rc)
            # clean old files
            i_fasta = move(i_fasta, tmp/(i_fasta.with_name(
                i_fasta.stem+'-noRC.fasta')).name)
            i_gb = move(i_gb, tmp/(i_gb.with_name(i_gb.stem+'-noRC.gb')).name)
            rc_fasta = move(rc_fasta, rc_fasta.with_suffix(''))
            r_rc_gb, r_rc_fasta = rotate_seq(rc_fasta, tmp=tmp)
            if r_rc_gb is None:
                continue
            rc_fasta.unlink()
            r_rc_gb = move(r_rc_gb, output/r_rc_gb.with_name(
                r_rc_gb.stem+'_RC.gb').name)
            r_rc_fasta = move(r_rc_fasta, output/r_rc_fasta.with_name(
                r_rc_fasta.stem+'_RC.fasta').name)
            new_compare_result = compare(r_rc_fasta, r_ref_fasta,
                                         arg.perc_identity)
            pdf = draw(r_ref_gb, r_rc_gb, new_compare_result)
            pdf = move(pdf, output/pdf.name)
            divided[i]['fasta'] = r_rc_fasta
            new_regions = get_regions(r_rc_gb)
            for _ in new_regions:
                divided[i][_] = len(new_regions[_])
            # validate again
            count_2, to_rc_2, *_ = validate_regions(
                option_len, new_regions, new_compare_result, arg.perc_identity)
            if to_rc_2 is None:
                success = True
        else:
            i_fasta = move(i_fasta, i_fasta.with_suffix('.fasta'))
            success = True
        divided[i]['success'] = success

    for i in divided:
        if divided[i]['success']:
            v_file = divided[i]['fasta']
            validated.append(v_file)
    if len(validated) != 0:
        log.info('Validated sequences:')
        for i in validated:
            log.info(f'\t{i.name}')
    output_info = output / f'{output.name}-results.csv'
    output_info_exist = output_info.exists()
    with open(output_info, 'a', encoding='utf-8') as out:
        if not output_info_exist:
            out.write('fasta,Success,Seed,Length,LSC,IRa,SSC,IRb,'
                      'Missing,Incomplete,RC_region,'
                      'Reference,Ref_length,r_LSC,r_IRa,r_SSC,r_IRb\n'
                      )
        for record in divided:
            # format is easier than f-string for dict
            simple = divided[record]
            # add seed info
            simple['seed'] = str(arg.seed)
            simple['fasta'] = simple['fasta'].name
            out.write('{fasta},{success},{seed},{length},{LSC},'
                      '{IRa},{SSC},{IRb},{missing},{incomplete},'
                      '{rc},'.format(**simple))
            out.write('{},{},{},{},{},{}\n'.format(
                r_ref_fasta.stem, ref_len, len(ref_regions['LSC']),
                len(ref_regions['IRa']), len(ref_regions['SSC']),
                len(ref_regions['IRb'])))
    log.info(f'Validation result was written into {output_info}')
    return validated, output_info


if __name__ == '__main__':
    validate_main()
