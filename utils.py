#!/usr/bin/python3

import argparse
import logging
from os import devnull, remove
from pathlib import Path
from subprocess import run
from tempfile import TemporaryDirectory
from time import sleep

from Bio import Entrez, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import reverse_complement as rc
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


# temporary directory
TMP = TemporaryDirectory()
NULL = open(devnull, 'w')
# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
TEMP_LOG = 'Temp.log'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO,
                    handlers=[logging.StreamHandler(),
                              logging.FileHandler(TEMP_LOG)])
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
log = logging.getLogger(__name__)


def down_ref(taxon, output):
    """
    Arg:
        taxon(str): given taxon name
        output(Path): output folder
    Return:
        ref(Path): gb file
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
        info = Entrez.read(Entrez.esummary(db='nuccore',
                                           webenv=handle['WebEnv'],
                                           query_key=handle['QueryKey']))
        accession = info[0]['Caption']
        count = int(handle['Count'])
        if count == 0:
            continue
        ref = output / f'{taxon}_{accession}.gb'
        content = Entrez.efetch(db='nuccore', webenv=handle['WebEnv'],
                                query_key=handle['QueryKey'], rettype='gb',
                                retmode='text', retmax=1)
        with open(ref, 'w') as out:
            out.write(content.read())
        return ref


def get_full_taxon(taxon):
    """
    Get full lineage of given taxon
    Return lineage list, only contains Kingdom, Phylum, Class, Order, Family,
    Genus, Species.
    Arg:
        taxon(str): given taxon name, could be common name, use quotation mark
        if necessary
    Return:
        lineage(list): lineage list, from higher to lower rank
    """
    split = taxon.split(' ')
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
    target = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus',
              'species']
    lineage = []
    for i in target:
        lineage.append(full_lineage.get(i, ''))
    if ' ' in lineage[-1]:
        lineage[-1] = lineage[-1].split(' ')[-1]
    # species name contains genus
    if lineage[-1] != '' and lineage[-2] != '':
        lineage[-1] = f'"{lineage[-2]} {lineage[-1]}"'
    return lineage


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
        exit(-1)
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


def repeat(filename, fmt):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.
    Assume the strand of sequence is in correct direction.
    Args:
        filename(Path): file to be repeat
        fmt(str): 'gb' or 'fasta'
    Return:
        repeated(Path): repeated file with origin format
        repeat_fasta(Path): fasta file for BLAST
    """
    new_file = filename.with_suffix('.repeat')
    new_fasta = filename.with_suffix('.repeat_fasta')
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
        seq(SeqRecord): gb record
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


def rotate_seq(filename, min_IR=1000):
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
    log.info(f'Rotate {filename}...')
    with open(filename, 'r') as _:
        start = _.readline()
        if start.startswith('>'):
            fmt = 'fasta'
        else:
            fmt = 'gb'
    # if fmt == 'fasta':
        # fasta = Path(filename)
        # gb = fasta.with_suffix('.gb')
        # SeqIO.convert(fasta, 'fasta', gb, 'gb', alphabet=IUPAC.ambiguous_dna)
    if fmt == 'gb':
        # gb = Path(filename)
        fasta = gb.with_suffix('.fasta')
        SeqIO.convert(gb, 'gb', fasta, 'fasta')
    else:
        fasta = filename
    new_fasta = fasta.with_suffix('.rotate')
    new_gb = fasta.with_suffix('.new_gb')
    success = False

    origin_seq = list(SeqIO.parse(filename, fmt))
    assert len(origin_seq) == 1
    origin_seq = origin_seq[0]
    origin_len = len(origin_seq)
    repeat_fasta = repeat(fasta)
    repeat_seq = SeqIO.read(repeat_fasta, 'fasta')
    blast_result = blast(repeat_fasta, repeat_fasta)
    # only one record, loop just for for unpack
    for query in parse_blast_tab(blast_result):
        locations = set()
        # use fasta's description
        name = ''
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
                log.critical(f'\tToo much difference. Reject.')
                continue
        if len(seq_IRa) != len(seq_IRb):
            log.warning(f'\tIRa ({len(seq_IRa)}) and IRb ({len(seq_IRb)}) do '
                        'not have same length.')
            if abs(len(seq_IRa) - len(seq_IRb)) > ambiguous_base_n:
                log.critical(f'\tToo much difference. Reject.')
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
        assert str(seq_SSC.seq) == str(
            new_seq.features[2].extract(new_seq).seq)
        success = True
        log.info(f'\tRegions of {name}:')
        log.info(f'\t\tLSC {new_seq.features[0].location}')
        log.info(f'\t\tIRa {new_seq.features[1].location}')
        log.info(f'\t\tSSC {new_seq.features[2].location}')
        log.info(f'\t\tIRb {new_seq.features[3].location}')
        with open(new_fasta, 'w') as out:
            out.write(f'>{name}\n{new_seq.seq}\n')
    with open(new_gb, 'w') as out3:
        raw_gb.features.extend(features)
        SeqIO.write(raw_gb, out3, 'gb')


    remove(blast_result)
    remove(repeat_fasta)
    if not success:
        log.critical(f'Failed to rotate {filename}.')
        raise SystemExit
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
    choices = ('LSC', 'IRa', 'SSC', 'IRb', 'whole')
    if choice not in choices:
        raise ValueError(f'Region must be in {choices}.')
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
    new_name = '_r_' + raw.name
    new_file = Path(gb.parent, '_r_' + gb.stem + '.fasta')
    with open(new_file, 'w') as out:
        out.write(f'>{new_name}\n')
        out.write(f'{new_seq}\n')
    # hotfix

    n_gb, n_fasta = rotate_seq(new_file)
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
    Rotate genbank or fasta file.
    Only handle the first record in file.
    Also contains other utilities.
    """
    arg = parse_args()
    new_gb, new_fasta = rotate_seq(arg.filename, arg.min_ir)
    return


if __name__ == '__main__':
    main()
