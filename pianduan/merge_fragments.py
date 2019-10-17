#!/usr/bin/python3

from Bio import SeqIO
from utils import blast, parse_blast_tab
from pathlib import Path
from sys import argv
from random import shuffle


def link(contigs):
    """
    Reorder contigs by overlap.
    Assume there is only one path to link.
    Args:
        contigs(list(SeqRecord)): contigs
    Return:
        link_info(list(blast_result)): link info of contigs
    """
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
            if qend != origin_len:
                # print('qend !=origin_len', hit)
                if origin_len - qend > ambiguous_base_n:
                    continue
                else:
                    print('to be continue')
                    continue
            if sstrand == 'plus' and sstart != 1:
                # consider ambiguous base or not?
                if sstart > ambiguous_base_n:
                    continue
                else:
                    print('to be continue')
                    continue
            if sstrand == 'minus':
                print('minus', hit)
                continue
            overlap.append(hit)
    print(*overlap, sep='\n')
    link_info = []
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
            # [None, None] as head/tail
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
            link_info.append(scaffold)
            try:
                scaffold = overlap_dict.popitem()[0]
            except KeyError:
                break
    # remove [None, None]
    link_info = [i[1:-1] for i in link_info]
    return link_info


def merge_contigs(contigs, link_info):
    """
    Use overlap information of contigs to merge them.
    Arg:
        contigs(list(SeqRecord)): contigs
        link_info(list(blast_result)): link info of contigs
    Return:
        merged(list(SeqRecord)): list of merged sequences
    """
    contigs_d = {i.id: i for i in contigs}
    merged = []
    for links in link_info:
        seq = contigs_d[links[0][0]]
        for link in links[:-1]:
            down = contigs_d[link[1]]
            seq += down[link[9]:]
        # tail do not have link after itself
        tail_seq = contigs_d[links[-1][1]]
        seq += tail_seq[links[-1][9]:]
        seq.id = f'Merged_sequence {len(seq)}bp'
        merged.append(seq)
    return merged


def main():
    # if contigs cannot be merged, they may be useless
    contigs = []
    for f in argv[1:]:
        fasta = Path(f)
        for idx, record in enumerate(SeqIO.parse(fasta, 'fasta')):
            record.id = f'{fasta.stem}-{idx}'
            record.description = ''
            contigs.append(record)
    # print([i.id for i in contigs])
    shuffle(contigs)
    link_info = link(contigs)
    merged = merge_contigs(contigs, link_info)
    SeqIO.write(merged, Path(argv[1]).with_suffix('.merge'), 'fasta')
    # link_info = link(contigs)
    # merged = merge_contigs(contigs, link_info)
    # print(merged)


if __name__ == '__main__':
    main()
