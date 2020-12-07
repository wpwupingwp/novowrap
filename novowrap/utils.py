#!/usr/bin/python3

from pathlib import Path
from shutil import unpack_archive
from subprocess import DEVNULL, run
from threading import Thread
from time import sleep
from urllib.request import urlopen
from zipfile import ZipFile
import gzip
import logging
import platform

from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import reverse_complement as rc
from Bio.SeqFeature import SeqFeature, FeatureLocation


# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
log = logging.getLogger('novowrap')
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
log = logging.getLogger('novowrap')


def accessible(name: Path, type_: str):
    """
    Check given path is accessible or not.
    Given path does not exist.
    Args:
        name(Path): folder or file, absolute path
        type_(str): 'folder' or 'file'
    Return:
        ok(bool): accessible or not
    """
    p = Path(name)
    if type_ == 'folder':
        try:
            p.mkdir()
            p.rmdir()
            ok = True
        except PermissionError:
            ok = False
    elif type_ == 'file':
        try:
            p.touch()
            p.unlink()
            ok = True
        except PermissionError:
            ok = False
    else:
        log.critical(f'Illegal type: {type_}')
        ok = False
    return ok


def move(source: Path, dest, copy=False):
    """
    Move source to dest and return dest.
    If set "copy", copy source to dest instead of move.
    Because Path.rename could not move file across different filesystem or
    drive, have to use copy and delete to implement "move".
    Warning:
        This function does not check whether dest exists or not.
    Args:
        source(Path): old path
        dest(Path or str): new path
        copy(bool): copy or move
    Return:
        dest(Path): new path
    """
    source = Path(source).absolute()
    dest = Path(dest).absolute()
    # avoid useless copy
    # Path.samefile may raise FileNotFoundError
    if source == dest:
        log.debug(f'{source} and {dest} are same.')
    else:
        # read_bytes/write_bytes includes open, read/write and close steps
        dest.write_bytes(source.read_bytes())
        if not copy:
            source.unlink()
    return dest


def get_full_taxon(taxon: str):
    """
    Get full lineage of given taxon, return lineage list.
    Not only contains Kingdom, Phylum, Class, Order, Family, Genus, Species.
    Arg:
        taxon(str): given taxon name, could be common name
    Return:
        ok(bool): query ok or not
        lineage(list): lineage list, from lower rank to higher, empty for fail
    """
    # make copy
    name = taxon.capitalize()
    try:
        search = Entrez.read(Entrez.esearch(db='taxonomy', term=f'"{name}"'))
    except Exception:
        return False, []
    if search['Count'] == '0':
        if ' ' not in name:
            return False, []
        if ' ' in name:
            name = name.split(' ')[0]
            sleep(0.5)
            try:
                search = Entrez.read(Entrez.esearch(db='taxonomy',
                                                    term=f'"{name}"'))
            except Exception:
                return [], False
            if search['Count'] == '0':
                return False, []
    taxon_id = search['IdList'][0]
    try:
        record = Entrez.read(Entrez.efetch(db='taxonomy', id=taxon_id))[0]
    except Exception:
        return False, []
    # names = [i['ScientificName'] for i in record['LineageEx']]
    full_lineage = [(i['Rank'], i['ScientificName']) for i in
                    record['LineageEx']]
    full_lineage.append((record['Rank'], record['ScientificName']))
    return True, reversed(full_lineage)


def get_ref(taxon: str, out: Path, tmp=None, mt_mode=False,
            simple_validate=False):
    """
    Get reference gb file.
    Only one record will be retrieved.
    Arg:
        taxon(str): given taxon name
        out(Path): output folder
        tmp(None or Path): temp folder
    Return:
        ref(Path or None): gb file, None for fail
        ref_taxon(str or None): taxon of reference's, may not be same with
        given taxon
    """
    log.info(f'Try to get reference of {taxon} from NCBI Genbank.')
    Entrez.email = 'guest@example.org'
    # handle quotation mark
    taxon = taxon.strip('"')
    taxon = taxon.strip("'")
    try:
        Entrez.read(Entrez.esummary(db='taxonomy', id='9606'))
    except Exception:
        log.critical('Failed to query on NCBI Entrez server, please check '
                     'your Internet connection.')
        return None, None
    if tmp is None:
        tmp = out
    get_taxon_ok = False
    lineage = []
    for try_ in range(3):
        get_taxon_ok, lineage = get_full_taxon(taxon)
        if get_taxon_ok:
            break
    if not get_taxon_ok:
        log.critical(f'Cannot find {taxon} in NCBI Taxonomy database.')
        log.warning("Please check the taxon name's spell.")
        return None, None
    for taxon in lineage:
        rank, taxon_name = taxon
        if taxon_name == '':
            continue
        if rank == 'order':
            log.critical('The taxonomy of the reference is not close-related.')
            log.critical('The result may be incorrect.')
        taxon_name = taxon_name.replace(' ', '_')
        # Entrez has limitation on query frequency (3 times per second)
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
        sleep(0.5)
        if mt_mode:
            query = (f'''{taxon_name}[Organism] AND refseq[filter] '''
                     f'''AND mitochondrion[filter]''')
        else:
            query = (f'''{taxon_name}[Organism] AND refseq[filter] '''
                     f'''AND (chloroplast[filter] OR plastid[filter])''')
        log.debug(f'Query from NCBI Genbank:\t{query}')
        # seems nuccore is more stable than taxonomy database
        handle = Entrez.read(Entrez.esearch(db='nuccore', term=query,
                                            usehistory='y'))
        count = int(handle['Count'])
        if count == 0:
            continue
        info = Entrez.read(Entrez.esummary(db='nuccore',
                                           webenv=handle['WebEnv'],
                                           query_key=handle['QueryKey']))
        accession = info[0]['Caption']
        ref = out / f'{taxon_name}_{accession}.gb'
        content = Entrez.efetch(db='nuccore', webenv=handle['WebEnv'],
                                query_key=handle['QueryKey'], rettype='gb',
                                retmode='text', retmax=1)
        with open(ref, 'w') as _:
            _.write(content.read())
        r_gb, r_fasta = rotate_seq(ref, tmp=tmp,
                                   simple_validate=simple_validate)
        if r_gb is not None:
            r_gb.unlink()
            r_fasta.unlink()
            log.info(f'Got {ref.name} as reference.')
            return ref, taxon_name
        else:
            continue
    return None, None


def blast(query: Path, target: Path, perc_identity=70):
    """
    Use simple BLAST with special output format.
    Args:
        query(Path): query filename
        target(Path): target filename
        perc_identity(float): perc_identity for BLAST
    Return:
        blast_out(Path): blast result filename
        blast_log(Path): blast log

    """
    fmt = ('qseqid sseqid sstrand qlen slen length pident gapopen qstart qend '
           'sstart send')
    query = query.absolute()
    target = target.absolute()
    blast_out = query.with_suffix('.blast')
    blast_log = query.with_suffix('.blast_log')
    # use blastn -subject instead of makeblastdb
    ok, blastn = get_blast()
    if not ok:
        log.critical('Failed to run BLAST.')
        return None, None
    b_log = open(blast_log, 'w')
    b_run = run(f'{blastn} -query {query} -subject {target} -outfmt "7 {fmt}" '
                f'-out {blast_out} -strand both -perc_identity '
                f'{perc_identity}',
                shell=True, stdout=b_log, stderr=b_log)
    b_log.close()
    # remove makeblastdb result
    if b_run.returncode != 0:
        log.critical('Cannot run BLAST.')
        return None, blast_log
    return blast_out, blast_log


def parse_blast_tab(filename: Path):
    """
    Parse BLAST result (tab format).
    Return [qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen,
    qstart, qend, sstart, send]
    Arg:
        filename(Path): blast result file
    Return:
        line(list): parsed result
    """
    query = []
    # empty for empty
    if filename is None:
        return []
    # blastn use system's default encoding
    with open(filename, 'r') as raw:
        for line in raw:
            if line.startswith('# BLAST'):
                if query:
                    yield query
                query.clear()
            elif line.startswith('#'):
                pass
            else:
                line = line.strip().split('\t')
                # last is str
                line[3:] = [int(float(i)) for i in line[3:]]
                query.append(line)


def repeat(filename, fmt):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.
    Args:
        filename(Path): file to be repeat
        fmt(str): 'gb' or 'fasta'
    Return:
        new_file(Path): repeated fasta
    """
    new_file = filename.with_suffix('.rpt')
    raw = SeqIO.read(filename, fmt)
    # assume given sequence is a whole chloroplast genome, no more or less
    new = raw + raw
    log.debug(f'Before {len(raw)}; Repeat {len(new)}')
    SeqIO.write(new, new_file, 'fasta')
    return new_file


def slice_gb(seq, location):
    """
    Biopython will skip incomplete annotation.
    Use this function to keep those fragments.
    Arg:
        seq(SeqRecord): sequence
        location(FeatureLocation): location
    Return:
        new_seq(SeqRecord): new SeqRecord
    """
    start = location.start
    end = location.stop
    f_before = None
    f_after = None
    for f in seq.features:
        if start in f:
            # slice may be inner of a feature
            f_before = SeqFeature(id=f.id, type='part',
                                  location=FeatureLocation(
                                      start-start,
                                      min(f.location.end, end)-start,
                                      f.location.strand),
                                  qualifiers=f.qualifiers)
            f_before.qualifiers['raw_location'] = str(f.location)
            break
    for f in seq.features:
        if end in f:
            f_after = SeqFeature(id=f.id, type='part',
                                 location=FeatureLocation(
                                     max(f.location.start, start)-start,
                                     end-start,
                                     f.location.strand),
                                 qualifiers=f.qualifiers)
            f_after.qualifiers['raw_location'] = str(f.location)
            break
    # both are not None
    if all([f_before, f_after]):
        if (f_before.location == f_after.location and
                f_before.qualifiers == f_after.qualifiers):
            # remove duplicated feature
            f_after = None
    new_seq = seq[location]
    for i in (f_before, f_after):
        if i is not None:
            i.qualifiers['sliced'] = str(location)
            i.qualifiers['reference'] = seq.id
            new_seq.features.append(i)
    return new_seq


def get_fmt(filename):
    """
    Detect file format.
    Support gz, fasta and gb.
    Args:
        filename(Path or str): filename
    Return:
        fmt(str): 'fasta', 'gz', 'gb' or 'txt' (others)
    """
    with gzip.open(filename) as _:
        try:
            _.readline()
            return 'gz'
        except OSError:
            pass
    with open(filename, 'r') as _:
        peek = _.readline()
        if peek.startswith('>'):
            fmt = 'fasta'
        elif peek.upper().startswith('LOCUS'):
            fmt = 'gb'
        else:
            fmt = 'txt'
    return fmt


def rotate_seq(filename, min_ir=1000, tmp=None, silence=True,
               simple_validate=False, retry=False):
    """
    Rotate genbank or fasta record, from LSC (trnH-psbA) to IRa, SSC, IRb.
    Input file should only contains one record.
    Arg:
        filename(Path): genbank or filename
        min_IR: minimum IR length
        tmp(Path or None): temporary folder
        silence(bool): print debug info or not
        simple_validate(bool): use abnormal rotate or not
        retry(bool): retry with less restrictions
    Return:
        new_gb(Path): gb file
        new_fasta(Path): fasta file
    """
    # normal rotate
    # allow 10 bases difference
    _allowed_difference = 10
    filename = Path(filename).absolute()
    if tmp is None:
        tmp = filename.parent
    if silence:
        log.setLevel(logging.CRITICAL)
    log.info(f'Rotate {filename}...')
    fmt = get_fmt(filename)
    # get origin seq
    origin_seq = list(SeqIO.parse(filename, fmt))
    assert len(origin_seq) == 1
    origin_seq = origin_seq[0]
    origin_seq.seq.alphabet = IUPAC.ambiguous_dna
    origin_len = len(origin_seq)
    # get repeat seq
    repeat_fasta = repeat(filename, fmt)
    repeat_seq = SeqIO.read(repeat_fasta, 'fasta')
    log.setLevel(logging.INFO)
    blast_result, blast_log = blast(repeat_fasta, repeat_fasta)
    log.setLevel(logging.CRITICAL)
    if blast_result is None:
        log.setLevel(logging.INFO)
        return None, None
    # clean tmp files immediately, in case of exceptions that break clean step
    repeat_fasta.unlink()
    new_fasta = filename.with_suffix('.rotate')
    if filename.suffix == '.gb':
        new_gb = filename.with_suffix('.gb.gb')
    else:
        new_gb = filename.with_suffix('.gb')
    if simple_validate:
        # do not rotate
        SeqIO.write(origin_seq, new_fasta, 'fasta')
        SeqIO.write(origin_seq, new_gb, 'gb')
        log.debug('simple rotate completed.')
        log.setLevel(logging.INFO)
        return new_gb, new_fasta
    success = False
    # only one record, loop just for for unpack
    for query in parse_blast_tab(blast_result):
        locations = set()
        # use fasta's description
        # for the second round
        ambiguous_base_n = len(str(origin_seq.seq).strip('ATCGatcg'))
        if retry:
            ambiguous_base_n = max(ambiguous_base_n, _allowed_difference)
            log.info(f'Allowing {ambiguous_base_n} bases difference this time.')
        # only use hit of IR to IR
        max_aln_n = 0
        # because of ambiguous base, percent of identity bases may not be 100%
        if ambiguous_base_n != 0:
            log.warning(f'\tFound {ambiguous_base_n} ambiguous bases.')
            p_ident_min = int((1-(ambiguous_base_n/origin_len))*100)
        else:
            p_ident_min = 100
        for hit in query:
            (qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen,
             qstart, qend, sstart, send) = hit
            # name = seq.name
            # only self
            if qseqid != sseqid:
                continue
            # skip too short match
            if length < min_ir:
                continue
            # too long, origin to repeat
            if length >= origin_len:
                continue
            # allow few gaps
            if gapopen != 0:
                if gapopen > ambiguous_base_n:
                    continue
            if pident < p_ident_min:
                continue
            # mismatch
            location = tuple(sorted([qstart, qend, sstart, send]))
            # hit across origin and repeat
            if location[-1] - location[0] > origin_len:
                continue
            # self to self or self to repeat self
            if len(set(location)) != 4:
                continue
            # filter short hit
            # the longest is IR
            if length < max_aln_n:
                continue
            else:
                max_aln_n = length
            locations.add(location)
        if not locations:
            continue
        locations = list(locations)
        locations.sort(key=lambda x: x[1]-x[0], reverse=True)
        locations.sort(key=lambda x: x[0])
        if len(locations) < 2:
            continue
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
        seq_LSC = slice_gb(repeat_seq, region_LSC)
        seq_IRa = slice_gb(repeat_seq, region_IRa)
        seq_SSC = slice_gb(repeat_seq, region_SSC)
        seq_IRb = slice_gb(repeat_seq, region_IRb)
        new_seq = seq_LSC + seq_IRa + seq_SSC + seq_IRb
        # to be continue
        new_seq.seq.alphabet = IUPAC.ambiguous_dna
        if origin_len != len(new_seq):
            log.warning(f'\tOld and new sequences do not have save length.')
            log.info(f'Old: {origin_len}\tNew: {len(new_seq)}')
            if abs(origin_len - len(new_seq)) > ambiguous_base_n:
                log.info(f'\tToo much difference. Reject.')
                continue
        if len(seq_IRa) != len(seq_IRb):
            log.warning(f'\tIRa ({len(seq_IRa)}) and IRb ({len(seq_IRb)}) do '
                        'not have same length.')
            if abs(len(seq_IRa) - len(seq_IRb)) > ambiguous_base_n:
                log.warning(f'\tToo much difference. Reject.')
                continue
        # output
        offset = -1
        features = []
        for f_name, f in zip(('large single copy (LSC)',
                              'inverted repeat A (IRa)',
                              'small single copy (SSC)',
                              'inverted repeat B (IRb)'),
                             (region_LSC, region_IRa, region_SSC, region_IRb)):
            length = f.stop - f.start
            features.append(SeqFeature(
                FeatureLocation(offset+1, length+offset+1),
                type='misc_feature',
                qualifiers={'note': f_name, 'software': 'rotate_seq'},
                strand=1))
            offset += length
        new_seq.features.extend(features)
        # print feature
        log.info(f'\tRegions of rotated record:')
        for f in features:
            log.info(f'\t{f.qualifiers["note"][-4:-1]}: {f.location.start} to '
                     f'{f.location.end} ({len(f)} bp)')
        SeqIO.write(new_seq, new_gb, 'gb')
        SeqIO.write(new_seq, new_fasta, 'fasta')
        success = True
    blast_result.unlink()
    blast_log.unlink()
    log.setLevel(logging.INFO)
    if not success:
        log.debug('Retry with less restrictions...')
        if retry:
            log.critical(f'Failed to rotate {filename}.')
            return None, None
        else:
            rotate_seq(filename, min_ir, tmp, silence, simple_validate,
                       retry=True)
    return new_gb, new_fasta


def get_regions(gb):
    """
    Arg:
        gb(Path): rotate_seq generated gb file, only contains one record
    Return:
        region({name: SeqFeature}): region location dict, default is ''
    """
    ref_region = dict().fromkeys(['LSC', 'SSC', 'IRa', 'IRb'], '')
    for feature in SeqIO.read(gb, 'gb').features:
        if (feature.type == 'misc_feature' and
                feature.qualifiers.get('software', ['', ])[0] == 'rotate_seq'):
            key = feature.qualifiers['note'][0][-4:-1]
            value = feature
            ref_region[key] = value
    return ref_region


def rc_regions(gb, choice='whole'):
    """
    Reverse and complement given region of sequence.
    Args:
        gb(Path or str): rotate_seq generated gb file
        choice(str): region to be processed, must be in 'LSC', 'IRa', 'SSC',
        'IRb', 'whole'.
    Return:
        new_file(Path): reverse-complemented fasta
    """
    # choices = ('LSC', 'IRa', 'SSC', 'IRb', 'whole')
    raw = SeqIO.read(gb, 'gb')
    data = {}
    new_seq = ''
    regions = get_regions(gb)
    for r in regions:
        data[r] = regions[r].extract(raw).seq
    if choice != 'whole':
        data[choice] = rc(regions[choice].extract(raw.seq))
        new_seq = data['LSC']
        for i in ['IRa', 'SSC', 'IRb']:
            new_seq += data[i]
    else:
        new_seq = rc(raw.seq)
    new_name = '_RC_' + raw.name
    new_file = gb.with_suffix('.rc.rc')
    with open(new_file, 'w') as _:
        _.write(f'>{new_name}\n')
        _.write(f'{new_seq}\n')
    return new_file


def test_cmd(program, option='-v'):
    """
    Test given program and option is ok to run or not.
    Args:
        program(Path or str): program path, could be relative path if it can
        be found in $PATH or %PATH%
        option(str): option for program, usually use "-v" to show version to
        test the program
    Return:
        success(bool): success or not
    """
    test = run(f'{program} {option}', shell=True, stdout=DEVNULL,
               stderr=DEVNULL)
    success = True if test.returncode == 0 else False
    return success


def get_third_party():
    """
    Get third_party folder.
    If do not exist, create it.
    If cannot access, report.
    Return:
        success(bool): ok or not
        third_party(Path): absolute path of third_party folder
    """
    third_party = Path().home().absolute() / '.novowrap'
    success = False
    if not third_party.exists():
        log.debug(f'Create folder {third_party}')
        try:
            third_party.mkdir()
        except Exception:
            log.critical(f'Failed to create {third_party}.'
                         'Please contact the administrator.')
            return success, third_party
    if not accessible(third_party/'test', 'file'):
        log.critical(f'Failed to access {third_party}.'
                     f'Please contact the administrator.')
        return success, third_party
    success = True
    return success, third_party


def get_perl(third_party=None):
    """
    Linux and Mac have perl already.
    For Windows user, this function help to get perl.exe .
    Only support x86_64 or amd64 machine.
    Args:
        third_party(Path or None): path for install
    Return:
        perl(str): perl location, empty for fail
    """
    url = ('http://strawberryperl.com/download/5.30.1.1/'
           'strawberry-perl-5.30.1.1-64bit-portable.zip')
    # have to test for threading
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return ''
    # only Windows need it
    home_perl = third_party / 'strawberry_perl' / 'perl' / 'bin' / 'perl.exe'
    # use run instead of find_executable because the later only check if exist
    # and ignore if could run
    perl = 'perl'
    if test_cmd(perl):
        return perl
    if test_cmd(home_perl):
        return str(home_perl)
    if platform.system() != 'Windows':
        log.critical('Cannot find Perl. Please follow the link '
                     'https://www.perl.org/get.html to install.')
        return ''
    else:
        log.warning('Cannot find Perl. Try to install.')
        if '64' not in platform.machine():
            log.critical(f'Unsupport machine {platform.machine()}')
            return ''
        try:
            # file is 148mb, 148mb/3000s=~50kb/s, consider it's ok for
            # most of users
            log.info('Download may be slow. Please consider manualy install '
                     'Perl according to https://www.perl.org/get.html ')
            down = urlopen(url, timeout=3000)
        except Exception:
            log.critical('Cannot download Perl. Please try to manually '
                         ' install.')
            return ''
        zip_file = third_party / 'strawberry-perl.zip'
        with open(zip_file, 'wb') as out:
            out.write(down.read())
        folder = third_party / 'strawberry_perl'
        try:
            with ZipFile(zip_file) as z:
                z.extractall(folder)
        except Exception:
            log.critical('The file is damaged.')
            log.critical('Please check your Internet connection.')
            return ''
        # fixed path in zip file, should not be wrong
        assert test_cmd(home_perl)
        return str(home_perl)


def get_novoplasty(third_party=None):
    """
    Ensure perl and novoplasty is available.
    Return novoplasty's path or None.
    Args:
        third_party(Path or None): path for install
    Return:
        novoplasty(Path or None): path of perl file, None for fail
    """
    url = 'https://github.com/ndierckx/NOVOPlasty/archive/NOVOPlasty3.7.2.zip'
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return None
    novoplasty = (third_party / 'NOVOPlasty-NOVOPlasty3.7.2' /
                  'NOVOPlasty3.7.2.pl')
    if novoplasty.exists():
        log.debug('Found NOVOPlasty in third_party folder.')
        return novoplasty
    log.critical('Cannot find NOVOPlasty, try to download.')
    # license warning
    log.info('\tThe program assumes that user ACCEPT the license of '
             'NOVOPlasty. ')
    log.info('\tIF NOT, please terminate the running.')
    log.info('\tSee https://raw.githubusercontent.com/'
             'ndierckx/NOVOPlasty/ for details.')
    try:
        # 126kb/10s=12.6kb/s, enough for test
        _ = urlopen('https://github.com', timeout=10)
    except Exception:
        log.critical('Cannot connect to github.com.')
        log.critical('Please check your Internet connection.')
        return None
    log.info('Due to connection speed, may need minutes. Please be patient.')
    try:
        # file is ~18mb, 18mb/360s=50kb/s, consider it's ok for most of users
        down = urlopen(url, timeout=360)
    except Exception:
        log.critical('Cannot download NOVOPlasty.')
        log.critical('Please manually download it from '
                     'https://github.com/ndierckx/NOVOPlasty')
        return None
    zip_file = third_party / 'NOVOPlasty3.7.2.zip'
    with open(zip_file, 'wb') as out:
        out.write(down.read())
    # novoplasty's files have illegal character ":" which cannot be extract by
    # shutil.unpack_archive in Windows, have to use zipfile
    # windows and linux both use "/"
    try:
        with ZipFile(zip_file, 'r') as z:
            z.extractall(third_party)
    except Exception:
        log.critical('The file is damaged.')
        log.critical('Please check your Internet connection.')
        return None
    log.info(f'Got {novoplasty.stem}.')
    return novoplasty


def get_blast(third_party=None):
    """
    Get BLAST location.
    If BLAST was found, assume makeblastdb is found, too.
    If not found, download it.
    Args:
        third_party(Path or None): path for install
    Return:
        ok(bool): success or not
        blast(str): blast path
    """
    if third_party is None:
        third_party_ok, third_party = get_third_party()
        if not third_party_ok:
            return third_party_ok, ''
    home_blast = third_party / 'ncbi-blast-2.10.0+' / 'bin' / 'blastn'
    # in Windows, ".exe" can be omitted
    # win_home_blast = home_blast.with_name('blastn.exe')
    ok = False
    # older than 2.8.1 is buggy
    url = ('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.10.0/'
           'ncbi-blast-2.10.0+')
    urls = {'Linux': url+'-x64-linux.tar.gz',
            'Darwin': url+'-x64-macosx.tar.gz',
            'Windows': url+'-x64-win64.tar.gz'}
    blast = 'blastn'
    if test_cmd(blast, '-version'):
        ok = True
        return ok, blast
    if test_cmd(home_blast, '-version'):
        ok = True
        return ok, str(home_blast)
    log.warning('Cannot find NCBI BLAST, try to install.')
    log.info('According to Internet speed, may be slow.')
    try:
        # 50kb/10s=5kb/s, enough for test
        _ = urlopen('https://www.ncbi.nlm.nih.gov', timeout=10)
    except Exception:
        log.critical('Cannot connect to NCBI.')
        log.critical('Please check your Internet connection.')
        return ok, ''
    try:
        # file is 86-222mb, 222mb/3600s=60kb/s, consider it's ok for users all
        # over the world
        down = urlopen(urls[platform.system()], timeout=3600)
    except Exception:
        log.critical('Cannot download BLAST.')
        log.critical('Please manually download it from'
                     'ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/')
        return ok, ''
    down_file = third_party / 'BLAST_2.10.0.tar.gz'
    with open(down_file, 'wb') as out:
        out.write(down.read())
    try:
        unpack_archive(down_file, third_party)
    except Exception:
        log.critical('The file is damaged.')
        log.critical('Please check your connection.')
        return ok, ''
    assert test_cmd(home_blast, '-version')
    ok = True
    return ok, str(home_blast)


def get_all_third_party():
    """
    Use three threads to speed up.
    """
    log.info('Try to locate or install all third-party software.')
    third_party_ok, third_party = get_third_party()
    if not third_party_ok:
        return -1
    perl = Thread(target=get_perl, args=(third_party,), daemon=True)
    novoplasty = Thread(target=get_novoplasty, args=(third_party,),
                        daemon=True)
    blast = Thread(target=get_blast, args=(third_party,), daemon=True)
    perl.start()
    novoplasty.start()
    blast.start()
    perl.join()
    novoplasty.join()
    blast.join()
    log.info('Finished.')


if __name__ == '__main__':
    get_all_third_party()
