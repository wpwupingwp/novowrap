#!/usr/bin/python3

from collections import defaultdict
from pathlib import Path
from sys import argv
from random import shuffle
import argparse

from Bio import SeqIO
from graphviz import Digraph

from utils import blast, parse_blast_tab


def parse_args(arg_list=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('input', nargs='+', help='input filenames')
    arg.add_argument('-o', '-out', dest='out', help='output filename')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def get_overlap(contigs):
    """
    Get overlap by BLAST.
    Args:
        contigs(list(SeqRecord)): contigs
    Return:
        link_info(list(blast_result)): link info of contigs
    """
    merged = Path('tmp.merged_fasta')
    SeqIO.write(contigs, merged, 'fasta')
    blast_result = blast(merged, merged)
    overlap = []
    for query, sequence in zip(parse_blast_tab(blast_result), contigs):
        if len(query) == 0:
            continue
        ambiguous_base_n = len(sequence.seq.strip('ATCGatcg'))
        # print('p ident min', p_ident_min)
        p_ident_min = int((1-(ambiguous_base_n/len(sequence)))*100)
        for hit in query:
            (qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen,
             qstart, qend, sstart, send) = hit
            if qseqid == sseqid:
                continue
            if pident < p_ident_min:
                continue
            # only allow overlapped seq
            # nested seqs should be omit
            if qlen == length or slen == length:
                continue
            # subject must be object's downstream
            if sstrand == 'plus':
                if (qend == qlen and sstart == 1):
                    pass
                else:
                    # consider ambiguous base or not?
                    continue
                    print('to be continue')
            if sstrand == 'minus':
                if (qstart == 1 and send == 1) or (
                        qend == qlen and sstart == slen):
                    pass
                else:
                    continue
            overlap.append(hit)
    # qseqid-sseqid: hit
    # ignore duplicate of plus-plus or minus-minus!
    #return overlap
    overlap_d1 = {tuple(sorted(i[:2])): i for i in overlap}
    raw = {tuple([i[0], i[1], i[2]]): i for i in overlap}
    left = {tuple([i[0], i[1], i[2]]): i for i in overlap_d1.values()}
    print('omit')
    for i in raw.keys()-left.keys():
        print(raw[i])
    print('left')
    print(*list(left.values()), sep='\n')
    print('all, left, left2', len(overlap), len(overlap_d1), len(left))
    return list(overlap_d1.values())


def remove_minus(overlap, contigs):
    """
    Use upstream/downstream information to find out contigs that should be
    reverse-complement (paired minus).
    Args:
        overlap(list(blast_result)): link info of contigs
        contigs(list(SeqRecord)): contigs
    Return:
        no_minus(list(SeqRecord)): contigs without minus
        minus_contig(list(SeqRecord)): minus contigs
    """
    a = []
    for i in contigs:
        i_rc = i.reverse_complement(id='_RC_'+i.id)
        a.extend([i, i_rc])
    return a, a
    plus = set()
    minus = set()
    minus_contig = []
    contigs_d = {i.id: i for i in contigs}
    for i in overlap:
        up, down, strand, *_ = i
        if strand == 'minus':
            print()
            minus.add(down)
            # minus.update([up, down])
        else:
            plus.update([up, down])
    to_rc = minus - plus
    print('torc', *to_rc)
    for i in to_rc:
        i_rc = contigs_d[i].reverse_complement(id='_RC_'+contigs_d[i].id)
        minus_contig.append(i_rc)
        contigs_d[i] = i_rc
    no_minus = list(contigs_d.values())
    # return contigs, contigs
    print('remove or not')
    return no_minus, minus_contig


def cut_circle(overlap):
    """
    Find circle, cut.
    """
    no_circle = []
    return no_circle


def clean_link(overlap):
    # seems useless
    """
    Ensure each upstream has only one downstream.
    Remove link that cause short circuit.
    Arg:
        overlap(list(blast_result)): link info of contigs
    Return:
        cleaned_link(list(blast_result)): clean link
    """
    # up_id: hit
    up_dict = defaultdict(set)
    # down_id: hit
    down_dict = defaultdict(set)
    raw = {(i[0], i[1]): i for i in overlap}
    bad_link = set()
    # seems defaultdict(list) cause "RuntimeError: dictionary changed size
    # during iteration"
    for i in overlap:
        up_dict[i[0]].add(i[1])
        down_dict[i[1]].add(i[0])
    for up, down in up_dict.items():
        if len(down) == 1:
            continue
        down_down = set()
        for d in down:
            if d in up_dict:
                down_down.update(up_dict[d])
        # one's downstreams should not overlap each other
        bad_link.update({(up, i) for i in down_down & down})
        good_down = down - down_down
        if len(good_down) == 1:
            continue
        for i in good_down:
            continue
            print(raw[(up, i)])
    cleaned_link = [raw[i] for i in raw if i not in bad_link]
    print('all, two up, bad, clean ')
    print(len(overlap), len(up_dict),  len(bad_link), len(cleaned_link))
    return cleaned_link


def get_link(contigs):
    """
    Reorder contigs by overlap.
    Assume there is only one path to link.
    Args:
        contigs(list(SeqRecord)): contigs
    Return:
        contigs_no_minus(list(SeqRecord)): contigs that remove minus
        link(list(blast_result)): link info of contigs
    """
    global dot
    dot = Digraph(engine='dot', node_attr={'shape': 'box'})
    overlap = get_overlap(contigs)
    contigs_no_minus, minus_contigs = remove_minus(overlap, contigs)
    overlap_no_minus = get_overlap(contigs_no_minus)
    # remove orphan minus
    up_dict = defaultdict(set)
    # down_id: hit
    down_dict = defaultdict(set)
    for i in overlap_no_minus:
        up_dict[i[0]].add(i[1])
        down_dict[i[1]].add(i[0])
    degree = {}
    for up in up_dict:
        print(up, up_dict[up])
        degree[up] = [0, len(up_dict[up])]
    for down in down_dict:
        if down in degree:
            degree[down] = [len(down_dict[down]), degree[down][1]]
        else:
            degree[down] = [len(down_dict[down]), 0]
    for i in degree:
        dot.node(i, xlabel=str(degree[i]))
    for i in overlap_no_minus:
        dot.node(i[0])
        dot.node(i[1])
        if i[2] == 'plus':
            dot.edge(i[0], i[1], color='#999999')
        else:
            dot.edge(i[0], i[1], color='#999999', style='dashed', dir='both')
    overlap_no_minus = [i for i in overlap_no_minus if i[2] != 'minus']
    # remove short circuit
    overlap_clean = clean_link(overlap_no_minus)
    # for i in overlap_clean:
    #    dot.edge(i[0], i[1], color='green')
    links = []
    scaffold = []
    # assume each seq only occurs once
    overlap_dict = {i[0]: i for i in overlap_clean}
    up_dict = {i[0]: i for i in overlap_clean}
    down_dict = {i[1]: i for i in overlap_clean}
    try:
        scaffold.append(overlap_dict.popitem()[1])
    except KeyError:
        raise Exception
    while True:
        try:
            up_name = scaffold[0][0]
            down_name = scaffold[-1][1]
        # name is None
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
                scaffold = [[None, None], *scaffold]
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
        else:
            scaffold.append([None, None])
        if scaffold[0][0] is None and scaffold[-1][0] is None:
            links.append(scaffold)
            try:
                scaffold = [overlap_dict.popitem()[1], ]
            except KeyError:
                break
    # remove [None, None]
    links = [i[1:-1] for i in links]
    for i in links:
        for j in i:
            dot.edge(j[0], j[1], color='blue')
    dot.render('graph')
    return contigs_no_minus, links


def merge_seq(contigs, links):
    """
    Use overlap information of contigs to merge them.
    Assume link_info only contains one kind of link route.
    Arg:
        contigs(list(SeqRecord)): contigs
        links(list(blast_result)): link info of contigs
    Return:
        merged(list(SeqRecord)): list of merged sequences
    """
    contigs_d = {i.id: i for i in contigs}
    merged = []
    for link in links:
        # qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen, qstart,
        # qend, sstart, send
        seq = contigs_d[link[0][0]]
        for node in link[:-1]:
            down = contigs_d[node[1]]
            seq += down[node[11]:]
        # tail do not have link after itself
        tail_seq = contigs_d[link[-1][1]]
        seq += tail_seq[link[-1][11]:]
        seq.id = f'Merged_sequence {len(seq)}bp'
        merged.append(seq)
    return merged


def merge_contigs(arg_str=None):
    if arg_str is None:
        arg = parse_args()
    else:
        arg = parse_args(arg_str.split(' '))
    if arg.out is None:
        arg.out = Path(arg.input[0]).with_suffix('.merge')
    contigs = []
    for f in arg.input:
        fasta = Path(f)
        for idx, record in enumerate(SeqIO.parse(fasta, 'fasta')):
            record.id = f'{fasta.stem}-{idx}'
            record.description = ''
            contigs.append(record)
    shuffle(contigs)
    contigs_no_minus, links = get_link(contigs)
    merged = merge_seq(contigs_no_minus, links)
    SeqIO.write(merged, Path(argv[1]).with_suffix('.merge'), 'fasta')


if __name__ == '__main__':
    merge_contigs()
