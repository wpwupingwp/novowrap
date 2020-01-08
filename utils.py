#!/usr/bin/python3

import argparse
import gzip
import logging
from os import devnull
from pathlib import Path
from subprocess import run
from time import sleep

from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import reverse_complement as rc
from Bio.SeqFeature import SeqFeature, FeatureLocation


NULL = open(devnull, 'w')
log = logging.getLogger('novowrap')


def move(source, dest, copy=False):
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
    if source.samefile(dest):
        pass
    else:
        # read_bytes/write_bytes includes open, read/write and close steps
        dest.write_bytes(source.read_bytes())
        if not copy:
            source.unlink()
    return Path(dest)


def get_full_taxon(taxon):
    """
    Get full lineage of given taxon, return lineage list.
    Not only contains Kingdom, Phylum, Class, Order, Family, Genus, Species.
    Arg:
        taxon(str): given taxon name, could be common name, use quotation mark
        if necessary
    Return:
        lineage(list): lineage list, from lower rank to higher
    """
    split = taxon.split(' ')
    # "Genus species var. blabla"
    if len(split) >= 2 and split[0][0].isupper() and split[0][1].islower():
        name = ' '.join(split[0:2])
    else:
        name = split[0]
    Entrez.email = 'guest@example.org'
    search = Entrez.read(Entrez.esearch(db='taxonomy', term=f'"{name}"'))
    if search['Count'] == '0':
        if ' ' not in name:
            log.critical(f'Cannot find {name} in NCBI Taxonomy database.')
            log.warning("Please check the taxon name's spell or the Internet "
                        "connection or NCBI server's status.")
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
    # names = [i['ScientificName'] for i in record['LineageEx']]
    full_lineage = [(i['Rank'], i['ScientificName']) for i in
                    record['LineageEx']]
    full_lineage.append((record['Rank'], record['ScientificName']))
    return reversed(full_lineage)


def get_ref(taxon, out):
    """
    Get reference gb file.
    Only one record will be retrieved.
    Arg:
        taxon(str): given taxon name
        out(Path): output folder
    Return:
        ref(Path): gb file
        ref_taxon(str): taxon of reference's, may not be same with given taxon
    """
    log.info(f'Try to get reference of {taxon} from NCBI Genbank.')
    lineage = get_full_taxon(taxon)
    if lineage is None:
        return None, None
    for taxon in lineage:
        rank, taxon_name = taxon
        if taxon_name == '':
            continue
        if rank == 'order':
            log.critical('The taxonomy of the reference is not close-related.')
            log.critical('The result may be incorrect.')
        taxon_name = taxon_name.replace(' ', '_')
        # Entrez has limitation on query frenquency (3 times per second)
        # https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
        sleep(0.5)
        query = (f'''{taxon_name}[Organism] AND refseq[filter] '''
                 f'''AND (chloroplast[filter] OR plastid[filter])''')
        log.debug(f'Query from NCBI Genbank:\t{query}')
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
        with open(ref, 'w', encoding='utf-8') as out:
            out.write(content.read())
        r_gb, r_fasta = rotate_seq(ref)
        if r_gb is None:
            continue
        else:
            r_gb.unlink()
            r_fasta.unlink()
            log.info(f'Got {ref.name} as reference.')
            return ref, taxon_name
    return None, None


def blast(query, target, perc_identity=70):
    """
    Use simple BLAST with special output format.
    Args:
        query(Path): query filename
        target(Path): target filename
        perc_identity(float): perc_identity for BLAST
    Return:
        blast_out(Path): blast result filename
    """
    fmt = ('qseqid sseqid sstrand qlen slen length pident gapopen qstart qend '
           'sstart send')
    blast_out = query.with_suffix('.blast')
    # use blastn -subject instead of makeblastdb
    b_run = run(f'blastn -query {query} -subject {target} -outfmt "7 {fmt}" '
                f'-out {blast_out} -strand both -perc_identity '
                f'{perc_identity}',
                shell=True, stdout=NULL, stderr=NULL)
    # remove makeblastdb result
    if b_run.returncode != 0:
        log.critical('Cannot run BLAST.')
        return None
    return blast_out


def parse_blast_tab(filename):
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
                # last is str
                line[3:] = [int(float(i)) for i in line[3:]]
                query.append(line)


def _repeat(filename, fmt):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.

    The repeated gb record contains two "source"-type features, which may not
    be parsed correctly by other functions, hence make this function private.
    Args:
        filename(Path): file to be repeat
        fmt(str): 'gb' or 'fasta'
    Return:
        new_file(Path): repeated file with origin format
        new_fasta(Path): fasta file for BLAST
    """
    new_file = filename.with_suffix('.rpt')
    new_fasta = filename.with_suffix('.rpt_fasta')
    raw = SeqIO.read(filename, fmt)
    # assume given sequence is a whole chloroplast genome, no more or less
    new = raw + raw
    SeqIO.write(new, new_file, fmt)
    if fmt != 'fasta':
        SeqIO.write(new, new_fasta, 'fasta')
    else:
        new_fasta = new_file
    return new_file, new_fasta


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
    with open(filename, 'r', encoding='utf-8') as _:
        peek = _.readline()
        if peek.startswith('>'):
            fmt = 'fasta'
        elif peek.upper().startswith('LOCUS'):
            fmt = 'gb'
        else:
            fmt = 'txt'
    return fmt


def rotate_seq(filename, min_ir=1000, silence=True):
    """
    Rotate genbank or fasta record, from LSC (trnH-psbA) to IRa, SSC, IRb.
    Input file should only contains one record.
    Arg:
        filename(Path or str): genbank or filename
        min_IR: minimum IR length
    Return:
        success(bool): success or not
    """
    if silence:
        log.setLevel(logging.CRITICAL)
    log.info(f'Rotate {filename}...')
    filename = Path(filename)
    fmt = get_fmt(filename)
    # get origin seq
    origin_seq = list(SeqIO.parse(filename, fmt))
    assert len(origin_seq) == 1
    origin_seq = origin_seq[0]
    origin_len = len(origin_seq)
    # get repeat seq
    repeat_record, repeat_fasta = _repeat(filename, fmt)
    repeat_seq = SeqIO.read(repeat_record, fmt)
    blast_result = blast(repeat_fasta, repeat_fasta)
    new_fasta = filename.with_suffix('.rotate')
    if filename.suffix == '.gb':
        new_gb = filename.with_suffix('.gb.gb')
    else:
        new_gb = filename.with_suffix('.gb')
    success = False
    # only one record, loop just for for unpack
    for query in parse_blast_tab(blast_result):
        locations = set()
        # use fasta's description
        ambiguous_base_n = len(str(origin_seq).strip('ATCGatcg'))
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
    repeat_record.unlink()
    if repeat_fasta.exists():
        repeat_fasta.unlink()
    log.setLevel(logging.INFO)
    if not success:
        log.critical(f'Failed to rotate {filename}.')
        return None, None
    return new_gb, new_fasta


def get_regions(gb):
    """
    Arg:
        gb(Path): rotate_seq generated gb file, only contains one record
    Return:
        region({name: SeqFeature}): region location dict
    """
    ref_region = {}
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
    with open(new_file, 'w', encoding='utf-8') as out:
        out.write(f'>{new_name}\n')
        out.write(f'{new_seq}\n')
    return new_file


def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=main.__doc__)
    arg.add_argument('filename', help='Genbank or fasta file')
    arg.add_argument('-min_ir', type=int, default=1000,
                     help='minimum length of IR region')
    return arg.parse_args()


def main():
    """
    Utilities for novowrap.
    If run directly, rotate input file.
    """
    arg = parse_args()
    new_gb, new_fasta = rotate_seq(arg.filename, arg.min_ir)
    return


if __name__ == '__main__':
    main()
