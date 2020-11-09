#!/usr/bin/python3

from os import environ
from pathlib import Path
from random import randint
import argparse
import logging

from Bio import SeqIO
from matplotlib import use as mpl_use
if environ.get('DISPLAY', '') == '':
    mpl_use('Agg')
from matplotlib import pyplot as plt
import numpy as np

from novowrap import utils


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
    arg.add_argument('-input', required=True, help='input filename')
    arg.add_argument('-ref', help='reference gb')
    arg.add_argument('-taxon', nargs='*', help='Taxonomy name')
    options = arg.add_argument_group('Option')
    options.add_argument('-simple_validate', action='store_true',
                         help='for plastids with abnormal structure')
    options.add_argument('-mt', dest='mt_mode', action='store_true',
                         help='for mitochondria (EXPERIMENTAL)')
    options.add_argument('-perc_identity', type=float, default=0.7,
                         help='minimum percentage of identity of BLAST, 0-100')
    options.add_argument('-len_diff', type=float, default=0.2,
                         help='maximum percentage of length difference of '
                         'query to reference, 0-100')
    options.add_argument('-seed',
                         help='seed used in assembly, only for caller')
    options.add_argument('-out', help='output folder')
    options.add_argument('-debug', action='store_true',
                         help='print debug info')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def init_arg(arg):
    """
    Initialize working folder with arg.
    Args:
        arg(NameSpace): arg generated from parse_args
    Return:
        success(bool): success or not
        arg(NameSpace): arg with output info
    """
    if arg.debug:
        logging.basicConfig(level=logging.DEBUG)
        try:
            import coloredlogs
            coloredlogs.install(level=logging.DEBUG, fmt=FMT, datefmt=DATEFMT)
        except ImportError:
            pass
    success = False
    arg.input = Path(arg.input).absolute()
    if not arg.input.exists():
        log.critical(f'Input file {arg.input} does not exist.')
        return success, arg
    if arg.ref is None and arg.taxon is None:
        log.warning('Nor reference either taxonomy was given.')
        return success, arg
    elif arg.ref is not None and arg.taxon is not None:
        log.critical('Cannot use "-taxon" and "-ref" at same time.')
        return success, arg
    if arg.ref is not None:
        arg.ref = Path(arg.ref).absolute()
    if arg.taxon is None:
        pass
    elif len(arg.taxon) > 1:
        # for "Genus species var. blabla", ignore subspecies words
        arg.taxon = ' '.join(arg.taxon[:2])
    else:
        arg.taxon = arg.taxon[0]
    if arg.out is None:
        arg.out = Path(arg.input.stem+'-out').absolute()
    else:
        arg.out = Path(arg.out).absolute()
    if __name__ != '__main__':
        # if called from novowrap, out exist and accessible already
        pass
    else:
        if arg.out.exists():
            log.warning(f'Output folder {arg.out.name} exists.')
            # give users one more chance
            new_name = arg.out.name + f'-{randint(0, 2020):04d}'
            arg.out = arg.out.with_name(new_name)
            if arg.out.exists() or not utils.accessible(arg.out, 'folder'):
                log.critical('Cannot create output folder.')
                return success, arg
            else:
                arg.out.mkdir()
                log.info(f'Use {arg.out.name} instead.')
        else:
            if not utils.accessible(arg.out, 'folder'):
                log.critical(f'Failed to access output folder {arg.out}.'
                             f'Please contact the administrator.')
                return success, arg
            arg.out.mkdir()
    arg.tmp = arg.out / 'Temp'
    # if called by novowrap, tmp is accessible already
    if not arg.tmp.exists():
        arg.tmp.mkdir()
    if arg.mt_mode:
        arg.simple_validate = True
    success = True
    return success, arg


def process_ref(arg):
    """
    Preprocess reference
    Prefer ref because assembly passed ref
    Returns:
        r_ref_gb(Path): rotated gb
        r_ref_fasta(Path): rotated fasta
        ref_len(int): length of the reference
    """
    r_ref_gb = None
    r_ref_fasta = None
    ref_len = 0
    if arg.ref is not None:
        log.info(f'Reference:\t{arg.ref}')
        fmt = utils.get_fmt(arg.ref)
        ref_gb = Path(arg.ref).absolute()
        ref_gb = utils.move(ref_gb, arg.tmp/ref_gb.name, copy=True)
        ref_records = list(SeqIO.parse(ref_gb, fmt))
        if len(ref_records) > 1:
            log.warning('Given reference contains more than one records, '
                        'only use the first.')
            # assume given reference is ok since user refuse to auto download
            # reference from Genbank
            SeqIO.write(ref_records[0], ref_gb, fmt)
    else:
        log.info(f'Taxonomy:\t{arg.taxon}')
        ref_gb, ref_taxon = utils.get_ref(arg.taxon, arg.tmp)
        if ref_gb is None:
            log.critical('Failed to get reference.')
            log.debug(f'{arg.input} {arg.ref} REF_NOT_FOUND\n')
            return r_ref_gb, r_ref_fasta, ref_len
        ref_gb = utils.move(ref_gb, arg.tmp/ref_gb.name)
        fmt = 'gb'
    log.info(f'Output:\t {arg.out}')
    ref_len = len(SeqIO.read(ref_gb, fmt))
    r_ref_gb, r_ref_fasta = utils.rotate_seq(
        ref_gb, tmp=arg.tmp, simple_validate=arg.simple_validate)
    if r_ref_gb is None:
        log.critical('Cannot process reference sequence.')
        log.critical('Please consider to use another reference.')
        log.debug(f'{arg.input} {arg.ref} REF_CANNOT_ROTATE\n')
        return r_ref_gb, r_ref_fasta, ref_len
    return r_ref_gb, r_ref_fasta, ref_len


def divide_records(fasta: Path, output: Path, ref_len: int,
                   tmp: Path, len_diff=0.1, simple_validate=False):
    """
    Make sure each file has only one record.
    Args:
        fasta(Path): fasta file
        output(Path): output folder
        ref_len(int): length of reference, to filter bad records
        tmp(Path): temp folder
        len_diff: maximum allowed length difference
        simple_validate(bool): skip normal rotate or not
    Returns:
        divided(dict): info of divided files
    """
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
            filename = output / f'{fasta.stem}_{idx+1}{fasta.suffix}'
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
            r_gb, r_fasta = utils.rotate_seq(
                filename, tmp=tmp, simple_validate=simple_validate)
            if r_gb is not None:
                divided[filename].update({'gb': r_gb, 'fasta': r_fasta,
                                          'length': record_len})
                utils.move(filename, tmp/filename.with_suffix('.raw').name)
            else:
                skip = 'structure_unusual'
        divided[filename]['skip'] = skip
        log.debug(f'{skip}')
    return divided


def compare_seq(query, reference, tmp, perc_identity):
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
    blast_result, blast_log = utils.blast(Path(query), reference,
                                          perc_identity*100)
    if blast_result is None:
        return None
    # only one record in file, loop is for unpack
    for query in utils.parse_blast_tab(blast_result):
        record = []
        for i in query:
            (qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen,
             qstart, qend, sstart, send) = i
            record.append([qstart, qend, sstart, send, sstrand, pident])
        results.append(record)
    utils.move(blast_result, tmp/blast_result.name)
    utils.move(blast_log, tmp/blast_log.name)
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


def draw(ref_gb: Path, seq_gb: Path, data: list):
    """
    Draw figure.
    Args:
        ref_gb(Path): reference genbank file
        seq_gb(Path): sequence genbank file
        data(list): BLAST result
    Return:
        pdf(Path): figure file
    """
    ref_regions = utils.get_regions(ref_gb)
    seq_regions = utils.get_regions(seq_gb)
    title = f'{seq_gb.stem} and {ref_gb.stem}'
    ignore_offset = len(ref_regions['IRa'])*2 + len(ref_regions['SSC'])
    plt.rcParams.update({'font.size': 20, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f'Validation of {title}', pad=30)
    plt.xlabel('Base')
    y_max = max(ref_regions['IRb'].location.end,
                seq_regions['IRb'].location.end)
    for key, value in ref_regions.items():
        plt.plot([value.location.start, value.location.end], [0.8, 0.8],
                 marker='+', label=key, linewidth=10)
        plt.text(value.location.start+len(value)/2, 0.78, f'{len(value)} bp',
                 fontsize=20, ha='center')
    for key, value in seq_regions.items():
        plt.text(value.location.start+len(value)/2, 0.96, f'{len(value)} bp',
                 fontsize=20, ha='center')
        if key != 'IRb':
            plt.axvline(x=value.location.end, ymin=0.15, ymax=0.85,
                        linestyle='dashdot', alpha=0.5)
    # no repeat legend
    plt.plot(0.5, 0.5, 'r-+', linewidth=5, label='Plus')
    plt.plot(0.5, 0.5, 'g-|', linewidth=5, label='Minus')
    plt.ylim([0.5, 1.1])
    plt.xlim([0, y_max])
    # usually fine
    if y_max > 20000:
        plt.xticks(list(range(0, y_max, 20000)))
    plt.yticks([0.65, 0.8, 0.95], labels=['Minus', 'Reference', 'Plus'])
    plt.legend(loc='lower left')
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


def simple_draw(ref_gb: Path, seq_gb: Path, data: list):
    """
    Draw figure of sequences without 4-parts structure.
    Args:
        ref_gb(Path): reference genbank file
        seq_gb(Path): sequence genbank file
        data(list): BLAST result
    Return:
        pdf(Path): figure file
    """
    ref_len = len(SeqIO.read(ref_gb, 'gb'))
    seq_len = len(SeqIO.read(seq_gb, 'gb'))
    title = f'{seq_gb.stem} and {ref_gb.stem}'
    ignore_offset = min(ref_len, seq_len)
    plt.rcParams.update({'font.size': 20, 'font.family': 'serif'})
    plt.figure(1, figsize=(30, 15))
    plt.title(f'Validation of {title}', pad=30)
    plt.xlabel('Base')
    y_max = max(ref_len, seq_len)
    plt.plot([0, ref_len], [0.8, 0.8], marker='+', label='Reference',
             linewidth=10)
    plt.text(ref_len/2, 0.78, f'{ref_len} bp', fontsize=20, ha='center')
    plt.text(ref_len/2, 0.96, f'{seq_len} bp', fontsize=20, ha='center')
    # no repeat legend
    plt.plot(0.5, 0.5, 'r-+', linewidth=5, label='Plus')
    plt.plot(0.5, 0.5, 'g-|', linewidth=5, label='Minus')
    plt.ylim([0.5, 1.1])
    plt.xlim([0, y_max])
    # usually fine
    if y_max > 20000:
        plt.xticks(list(range(0, y_max, 20000)))
    plt.yticks([0.65, 0.8, 0.95], labels=['Minus', 'Reference', 'Plus'])
    plt.legend(loc='lower left')
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


def validate_regions(length: int, regions: dict, compare: list,
                     perc_identity=0.7):
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


def write_output(divided, r_ref_fasta, ref_len, ref_regions, arg):
    """
    Write output info.
    Args:
        divided:
        r_ref_fasta:
        ref_len:
        ref_regions:
        arg:
    Returns:
        output_info(Path): csv file
    """
    output_info = arg.out / f'{arg.out.name}-results.csv'
    if not output_info.exists():
        with open(output_info, 'w') as csv_head:
            csv_head.write('Input,Success,Seed,Length,LSC,IRa,SSC,IRb,'
                           'Missing,Incomplete,RC_region,'
                           'Reference,Ref_length,r_LSC,r_IRa,r_SSC,r_IRb\n')
    with open(output_info, 'a') as out:
        for record in divided:
            # format is easier than f-string for dict
            simple = divided[record]
            # add seed info
            simple['seed'] = str(arg.seed)
            simple['fasta'] = simple['fasta'].stem
            out.write('{fasta},{success},{seed},{length},{LSC},'
                      '{IRa},{SSC},{IRb},{missing},{incomplete},'
                      '{rc},'.format(**simple))
            out.write('{},{},{},{},{},{}\n'.format(
                r_ref_fasta.stem, ref_len, len(ref_regions['LSC']),
                len(ref_regions['IRa']), len(ref_regions['SSC']),
                len(ref_regions['IRb'])))
    return output_info


def normal_validate(divided: dict, r_ref_gb: Path, r_ref_fasta: Path, arg):
    """
    For chloroplast genomes with 4-parts structure
    """
    for i in divided:
        success = False
        divided[i]['success'] = success
        if divided[i]['skip']:
            continue
        i_gb = divided[i]['gb']
        i_fasta = divided[i]['fasta']
        log.info(f'Analyze {i_fasta}.')
        option_regions = utils.get_regions(i_gb)
        # add regions info
        for _ in option_regions:
            divided[i][_] = len(option_regions[_])
        compare_result = compare_seq(i_fasta, r_ref_fasta, arg.tmp,
                                     arg.perc_identity)
        if compare_result is None:
            log.critical('Cannot run BLAST.')
            log.debug(f'{arg.input} {arg.ref} BLAST_FAIL\n')
            return None
        pdf = draw(r_ref_gb, i_gb, compare_result)
        utils.move(pdf, arg.out/pdf.name)
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
            rc_fasta = utils.rc_regions(i_gb, to_rc)
            # clean old files
            utils.move(i_fasta, arg.tmp/(i_fasta.with_name(
                i_fasta.stem+'-noRC.fasta')).name)
            utils.move(i_gb, arg.tmp/(i_gb.with_name(
                i_gb.stem+'-noRC.gb')).name)
            rc_fasta = utils.move(rc_fasta, rc_fasta.with_suffix(''))
            r_rc_gb, r_rc_fasta = utils.rotate_seq(
                rc_fasta, tmp=arg.tmp, simple_validate=arg.simple_validate)
            if r_rc_gb is None:
                continue
            rc_fasta.unlink()
            r_rc_gb = utils.move(r_rc_gb, arg.out/r_rc_gb.with_name(
                r_rc_gb.stem+'_RC.gb').name)
            r_rc_fasta = utils.move(r_rc_fasta, arg.out/r_rc_fasta.with_name(
                r_rc_fasta.stem+'_RC.fasta').name)
            new_compare_result = compare_seq(r_rc_fasta, r_ref_fasta, arg.tmp,
                                             arg.perc_identity)
            pdf = draw(r_ref_gb, r_rc_gb, new_compare_result)
            utils.move(pdf, arg.out/pdf.name)
            divided[i]['fasta'] = r_rc_fasta
            new_regions = utils.get_regions(r_rc_gb)
            for _ in new_regions:
                divided[i][_] = len(new_regions[_])
            # validate again
            count_2, to_rc_2, *_ = validate_regions(
                option_len, new_regions, new_compare_result, arg.perc_identity)
            if to_rc_2 is None:
                success = True
        else:
            utils.move(i_fasta, i_fasta.with_suffix('.fasta'))
            success = True
        divided[i]['success'] = success
        return divided


def simple_validate(divided: dict, r_ref_gb: Path, r_ref_fasta: Path, arg):
    """
    For chloroplast genomes WITHOUT 4-parts structure
    """
    for i in divided:
        success = False
        divided[i]['success'] = success
        if divided[i]['skip']:
            continue
        i_gb = divided[i]['gb']
        i_fasta = divided[i]['fasta']
        log.info(f'Analyze {i_fasta}.')
        option_regions = utils.get_regions(i_gb)
        # add regions info
        for _ in option_regions:
            divided[i][_] = len(option_regions[_])
        compare_result = compare_seq(i_fasta, r_ref_fasta, arg.tmp,
                                     arg.perc_identity)
        if compare_result is None:
            log.critical('Cannot run BLAST.')
            log.debug(f'{arg.input} {arg.ref} BLAST_FAIL\n')
            return None
        pdf = simple_draw(r_ref_gb, i_gb, compare_result)
        utils.move(pdf, arg.out/pdf.name)
        log.debug('Skip detecting reverse complement region.')
        utils.move(i_fasta, i_fasta.with_suffix('.fasta'))
        success = True
        divided[i]['success'] = success
        return divided


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
    init_ok, arg = init_arg(arg)
    if not init_ok:
        log.critical('Quit.')
        return validated, output_info
    log.info(f'Input:\t{arg.input}')
    r_ref_gb, r_ref_fasta, ref_len = process_ref(arg)
    if r_ref_gb is None:
        return validated, output_info

    ref_regions = utils.get_regions(r_ref_gb)
    divided = divide_records(arg.input, arg.out, ref_len, arg.tmp,
                             arg.len_diff, arg.simple_validate)
    if arg.simple_validate:
        divided = simple_validate(divided, r_ref_gb, r_ref_fasta, arg)
    else:
        divided = normal_validate(divided, r_ref_gb, r_ref_fasta, arg)
    if divided is None:
        return validated, output_info
    for i in divided:
        if divided[i]['success']:
            v_file = divided[i]['fasta']
            validated.append(v_file)
    if len(validated) != 0:
        log.info('Validated sequences:')
        for i in validated:
            log.info(f'\t{i.name}')
    output_info = write_output(divided, r_ref_fasta, ref_len, ref_regions, arg)
    log.info(f'Validation result was written into {output_info}')
    return validated, output_info


if __name__ == '__main__':
    validate_main()
