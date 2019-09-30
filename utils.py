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
    Args:
        source(Path): old path
        dest(Path or str): new path
        copy(bool): copy or move
    Return:
        dest(Path): new path
    """
    if not copy:
        source.rename(dest)
        return Path(dest)
    elif source.absolute() == dest.absolute():
        return Path(dest)
    else:
        with open(source, 'rb') as a, open(dest, 'wb') as b:
            b.write(a.read())
        return Path(dest)


def get_full_taxon(taxon):
    """
    Get full lineage of given taxon
    Return lineage list, only contains Kingdom, Phylum, Class, Order, Family,
    Genus, Species.
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
    target = ['species', 'genus', 'family', 'order', 'class', 'phylum',
              'kingdom']
    lineage = []
    for i in target:
        lineage.append([i, full_lineage.get(i, '')])
    # if ' ' in lineage['species']:
    #     lineage['species'] = lineage['species'].split(' ')[-1]
    # species name contains genus
    # if lineage[-1] != '' and lineage[-2] != '':
    #     lineage[-1] = f'"{lineage[-2]} {lineage[-1]}"'
    return lineage


def get_ref(taxon):
    """
    Get reference gb file.
    Only one record will be retrieved.
    Arg:
        taxon(str): given taxon name
    Return:
        ref(Path): gb file
        accession(str): accession number of the record
    """
    lineage = get_full_taxon(taxon)
    if lineage is None:
        log.warning(f'Cannot find {taxon} in NCBI taxonomy database.')
        log.warning("Please check the taxon name's spell or the Internet "
                    "connection or NCBI server's status.")
        return None
    for taxon in lineage:
        rank, taxon_name = taxon
        if taxon_name == '':
            continue
        if rank == 'class':
            log.critical('Cannot find close related reference.')
            log.critical('The result may be incorrect.')
        # use underscore to replace space
        if ' ' in taxon_name:
            taxon_name = taxon_name.strip('"')
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
        ref = Path(f'{taxon_name}_{accession}.gb')
        content = Entrez.efetch(db='nuccore', webenv=handle['WebEnv'],
                                query_key=handle['QueryKey'], rettype='gb',
                                retmode='text', retmax=1)
        with open(ref, 'w') as out:
            out.write(content.read())
        return ref


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
    FMT = ('qseqid sseqid sstrand length pident gapopen qstart qend '
           'sstart send')
    blast_out = query.with_suffix('.blast')
    # use blastn -subject instead of makeblastdb
    blast = run(f'blastn -query {query} -subject {target} -outfmt "7 {FMT}" '
                f'-out {blast_out} -strand both -perc_identity '
                f'{perc_identity}',
                shell=True, stdout=NULL, stderr=NULL)
    # remove makeblastdb result
    if blast.returncode != 0:
        log.critical('Cannot run BLAST.')
        return None
    return blast_out


def parse_blast_tab(filename):
    """
    Parse BLAST result (tab format).
    Return [qseqid, sseqid, sstrand, length, pident, gapopen,
    qstart, qend, sstart, send]
    Arg:
        filename(Path): blast result file
    Return:
        line(list): parsed result
    """
    query = []
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
    with open(filename, 'r') as _:
        peek = _.readline()
        if peek.startswith('>'):
            fmt = 'fasta'
        elif peek.upper().startswith('LOCUS'):
            fmt = 'gb'
        else:
            fmt = 'txt'
    return fmt


def rotate_seq(filename, min_IR=1000, silence=True):
    """
    Rotate genbank or fasta record, from LSC (trnH-psbA) to IRa, SSC, IRb.
    Input file should only contains one record.
    Arg:
        filename(Path or str): genbank or filename
        min_IR: minimum IR length
    Return:
        success(bool): success or not
    """
    # FMT = 'qseqid sseqid length pident gapopen qstart qend sstart
    # send'
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
            (qseqid, sseqid, sstrand, length, pident, gapopen, qstart, qend,
             sstart, send) = hit
            # name = seq.name
            # only self
            if qseqid != sseqid:
                continue
            # skip too short match
            if length < min_IR:
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
        n_gb(Path): gb file after reverse-complement and rotate
        n_fasta(Path): fasta file after reverse-complement and rotate
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
    new_file = gb.with_name(gb.stem+'_RC'+'.fasta')
    with open(new_file, 'w') as out:
        out.write(f'>{new_name}\n')
        out.write(f'{new_seq}\n')
    # hide rotate log
    n_gb, n_fasta = rotate_seq(new_file)
    new_file.unlink()
    return n_gb, n_fasta


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
