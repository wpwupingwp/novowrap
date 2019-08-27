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


def down_ref(taxon):
    """
    Arg:
        taxon(str): given taxon name
    Return:
        output_file(Path): fasta file
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
        count = int(handle['Count'])
        if count == 0:
            continue
        output_file = f'{taxon.replace(" ", "_")}.gb'
        content = Entrez.efetch(db='nuccore', webenv=handle['WebEnv'],
                                query_key=handle['QueryKey'], rettype='gb',
                                retmode='text', retmax=1)
        output_file.write_text(content.read())
        return output_file



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
    FMT = ('qseqid sseqid qseq sseq sstrand qlen pident gapopen qstart qend '
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
    Return [qseqid, sseqid, qseq, sseq, sstrand, length, pident, gapopen,
    qstart, qend, sstart, send, sstrand]
    Arg:
        filename(Path): blast result file
    Return:
        line(list): parsed result
    """
    query = []
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
                line[5:] = list(map(float, line[5:]))
                line[5:] = list(map(int, line[5:]))
                query.append(line)


def repeat(fasta):
    """
    Duplicate sequence to get full length of IR which may locate in start or
    end of sequences that BLAST cannot get whole length.
    Assume the strand of sequence is in correct direction.
    Args:
        fasta(Path): fasta filename
        taxon(Path): taxonomy of given fasta
    Return:
        new_fasta(Path): new_fasta's name
    """
    new_fasta = fasta.with_suffix('.new')
    new = []
    for i in SeqIO.parse(fasta, 'fasta'):
        i.seq = i.seq + i.seq
        new.append(i)
    SeqIO.write(new, new_fasta, 'fasta')
    return new_fasta


def slice_gb(seq, location):
    """
    Biopython will skip incomplete annotation.
    Use this function to keep those fragments.
    Arg:
        seq(SeqRecord): gb record
        location(slice): location
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
    Repeat length parameters for easily call rotate function instead of run
    whole program.
    Arg:
        filename(Path or str): genbank or filename
        min_IR: minimum IR length
    Return:
        success(bool): success or not
    """
    # FMT = 'qseqid sseqid qseq sseq length pident gapopen qstart qend sstart
    # send'
    with open(filename, 'r') as _:
        start = _.readline()
        if start.startswith('>'):
            fmt = 'fasta'
        else:
            fmt = 'gb'
    records = list(SeqIO.parse(filename, fmt))
    if len(records) > 1:
        log.warning(f'Found {len(records)} records')
        log.warning(f'Only handle the first as representative.')
        log.warning(f'Divide them if you want to rotate all.')
        filename = str(filename) + '.1'
        SeqIO.write(records[0], filename, fmt)
    if fmt == 'fasta':
        fasta = Path(filename)
        gb = fasta.with_suffix('.gb')
        SeqIO.convert(fasta, 'fasta', gb, 'gb', alphabet=IUPAC.ambiguous_dna)
    else:
        gb = Path(filename)
        fasta = gb.with_suffix('.fasta')
        SeqIO.convert(gb, 'gb', fasta, 'fasta')

    new_fasta = fasta.with_suffix('.rotate')
    new_regions = fasta.with_suffix('.regions')
    new_gb = fasta.with_suffix('.new.gb')
    success = False
    repeat_fasta = repeat(fasta)
    blast_result = blast(repeat_fasta, repeat_fasta)
    # blast and fasta use same order
    for query, seq, raw_gb in zip(parse_blast_tab(blast_result),
                                  SeqIO.parse(repeat_fasta, 'fasta'),
                                  SeqIO.parse(gb, 'gb')):
        log.info(f'Analyze {seq.name}.')
        locations = set()
        # use fasta's description
        name = ''
        original_seq_len = len(seq) // 2
        ambiguous_base_n = len(str(seq).strip('ATCGatcg')) // 2
        # only use hit of IR to IR
        max_aln_n = 0
        # because of ambiguous base, percent of identity bases may not be 100%
        if ambiguous_base_n != 0:
            log.warning(f'Found {ambiguous_base_n} ambiguous bases.')
            p_ident_min = int((1-(ambiguous_base_n/original_seq_len))*100)
        else:
            p_ident_min = 100
        for hit in query:
            (qseqid, sseqid, qseq, sseq, sstrand, qlen, pident, gapopen,
             qstart, qend, sstart, send) = hit
            name = seq.name
            # only self
            if qseqid != sseqid:
                continue
            # skip too short match
            if len(qseq) < min_IR:
                continue
            # origin to repeat
            if len(qseq) == len(seq):
                continue
            # allow few gaps
            if gapopen != 0:
                log.warning(f'Found {gapopen} gaps.')
                if gapopen > ambiguous_base_n:
                    log.critical('Too much gaps. Reject.')
                    continue
            # mismatch
            location = tuple(sorted([qstart, qend, sstart, send]))
            if pident < p_ident_min:
                continue
            # hit across origin and repeat
            if location[-1] - location[0] > original_seq_len:
                continue
            # self to self or self to repeat self
            if len(set(location)) != 4:
                continue
            # filter short hit
            if len(qseq) < max_aln_n:
                continue
            else:
                max_aln_n = len(qseq)
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
        seq_LSC = slice_gb(seq, region_LSC)
        seq_IRa = slice_gb(seq, region_IRa)
        seq_SSC = slice_gb(seq, region_SSC)
        seq_IRb = slice_gb(seq, region_IRb)
        new_seq = seq_LSC + seq_IRa + seq_SSC + seq_IRb
        new_seq.seq.alphabet = IUPAC.ambiguous_dna
        if len(seq)//2 != len(new_seq):
            log.warning(f'Old and new sequences do not have save length.')
            log.info(f'Old: {len(seq)//2}\tNew: {len(new_seq)}')
            if abs(len(seq)//2 - len(new_seq)) > ambiguous_base_n:
                log.critical(f'Too much difference. Reject.')
                continue
        if len(seq_IRa) != len(seq_IRb):
            log.warning(f'IRa ({len(seq_IRa)}) and IRb ({len(seq_IRb)}) do '
                        'not have same length! Reject.')
            if abs(len(seq_IRa) - len(seq_IRb)) > ambiguous_base_n:
                log.critical(f'Too much difference. Reject.')
                continue
        # output
        offset = -1
        features = []
        # for ogdraw
        for f_name, f in zip(('large single copy (LSC)',
                              'inverted repeat A (IRa)',
                              'small single copy (SSC)',
                              'inverted repeat B (IRb)'),
                             (region_LSC, region_IRa, region_SSC, region_IRb)):
            length = f.stop - f.start
            features.append(SeqFeature(
                FeatureLocation(offset+1, length+offset+1),
                type='misc_feature',
                qualifiers={'note': f_name, 'software': 'rotate_gb'},
                strand=1))
            offset += length
        # output
        new_seq.features.extend(features)
        assert str(seq_SSC.seq) == str(
            new_seq.features[2].extract(new_seq).seq)
        log.info(f'Rotated regions of {name}:')
        log.info(f'\tLSC {new_seq.features[0].location}')
        log.info(f'\tIRa {new_seq.features[1].location}')
        log.info(f'\tSSC {new_seq.features[2].location}')
        log.info(f'\tIRb {new_seq.features[3].location}')
        with open(new_fasta, 'w') as out:
            out.write(f'>{name}\n{new_seq.seq}\n')
        with open(new_regions, 'w') as o:
            o.write(f'>{name}-LSC\n{seq_LSC.seq}\n')
            o.write(f'>{name}-IRa\n{seq_IRa.seq}\n')
            o.write(f'>{name}-SSC\n{seq_SSC.seq}\n')
            o.write(f'>{name}-IRb\n{seq_IRb.seq}\n')
        with open(new_gb, 'w') as out3:
            raw_gb.features.extend(features)
            SeqIO.write(raw_gb, out3, 'gb')
        success = True

    remove(blast_result)
    remove(repeat_fasta)
    if not success:
        return None, None, None
    log.info(f'Rotated {gb.name} to uniform conformation {new_gb.name}.')
    return new_gb, new_fasta, new_regions


def rc(fasta, region=None, choice='whole'):
    """
    Reverse and complement given region of sequence.
    Args:
        fasta(Path or str): fasta file
        region(Path or str): fasta file contains regions
        choice(str): region to be processed, must be in 'LSC', 'IRa', 'SSC',
        'IRb', 'whole'.
    Return:
        rc(str): processed file
    """
    choices = ('LSC', 'IRa', 'SSC', 'IRb', 'whole')
    if choice not in choices:
        raise ValueError(f'Region must be in {choices}.')
    if region is None:
        # unrotated
        (r_gb, r_fasta, r_regions) = rotate_seq(fasta)
        region = r_regions
    raw = SeqIO.read(fasta, 'fasta')
    new_name = raw.name + '_rc'
    new_file = fasta + '.rc'
    data = {}
    for i in SeqIO.parse(r_regions, 'fasta'):
        name = i.id.split('-')[-1]
        data[name] = i
    if choice != 'whole':
        data[choice] = data[choice].reverse_complement(
            id=new_name+'_rc')
        new = data['LSC']
        for i in ['IRa', 'SSC', 'IRb']:
            new += data[i]
        SeqIO.write(new, new_file, 'fasta')
    else:
        SeqIO.write(raw.reverse_complement(id=new_name),
                    new_file, 'fasta')


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
    (new_gb, new_fasta, new_regions) = rotate_seq(arg.filename, arg.min_ir)
    return


if __name__ == '__main__':
    main()
