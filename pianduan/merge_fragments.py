#!/usr/bin/python3

from Bio import SeqIO
from utils import blast, parse_blast_tab
from pathlib import Path
from sys import argv
from random import shuffle


def lianjie(contigs):
    merged = Path('tmp.merged_fasta')
    SeqIO.write(contigs, merged, 'fasta')
    blast_result = blast(merged, merged)
    overlap = []
    for query, seq in zip(parse_blast_tab(blast_result), contigs):
        ambiguous_base_n = len(seq.seq.strip('ATCGatcg'))
        origin_len = len(seq)
        p_ident_min = int((1-(ambiguous_base_n/origin_len))*100)
        for hit in query:
            (qseqid, sseqid, sstrand, length, pident, gapopen, qstart, qend,
             sstart, send) = hit
            if qseqid == sseqid:
                continue
            if pident < p_ident_min:
                continue
            # only allow overlapped seq
            if qend != origin_len or sstart != 1:
                if sstart > ambiguous_base_n:
                    continue
                else:
                    print()
                    pass
            overlap.append(hit)
    print('*'*80)
    genome = []
    scaffold = []
    overlap_dict = {}
    # assume each seq only occurs once
    overlap_dict = {i[0]: i for i in overlap}
    up_dict = {i[0]: i for i in overlap}
    down_dict = {i[1]: i for i in overlap}
    scaffold.append(overlap_dict.popitem()[1])
    while True:
        try:
            up_name = scaffold[0][0]
            down_name = scaffold[-1][1]
        except TypeError:
            pass
        # if downstream not in overlap_dict:
        if up_name is None:
            pass
        elif up_name in down_dict:
            upstream = down_dict[up_name][0]
            if upstream in overlap_dict:
                scaffold = [overlap_dict.pop(upstream), *scaffold]
            else:
                print(f'{upstream} was used more than once!')
                pass
        else:
            scaffold = [[None, None], *scaffold]
        if down_name is None:
            pass
        elif down_name in up_dict:
            downstream = up_dict[down_name][0]
            if downstream in overlap_dict:
                scaffold.append(overlap_dict.pop(downstream))
        else:
            scaffold.append([None, None])
        if scaffold[0][0] is None and scaffold[-1][0] is None:
            genome.append(scaffold)
            try:
                scaffold = overlap_dict.popitem()[0]
            except KeyError:
                break
    print('-'*80)
    print('genome', *genome, sep='\n')
    return scaffold


def main():
    # if contigs cannot be merged, they may be useless
    contigs = []
    for f in argv[1:]:
        fasta = Path(f)
        print(fasta)
        for idx, record in enumerate(SeqIO.parse(fasta, 'fasta')):
            record.id = f'{fasta.stem}-{idx}'
            record.description = ''
            contigs.append(record)
    lianjie(contigs)
    shuffle(contigs)
    lianjie(contigs)


if __name__ == '__main__':
    main()
