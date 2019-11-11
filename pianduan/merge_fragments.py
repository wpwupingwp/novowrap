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
    overlap_d1 = {tuple(sorted(i[:2])): i for i in overlap}
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
            minus.add(down)
            # minus.update([up, down])
        else:
            plus.update([up, down])
    to_rc = minus - plus
    for i in to_rc:
        i_rc = contigs_d[i].reverse_complement(id='_RC_'+contigs_d[i].id)
        minus_contig.append(i_rc)
        contigs_d[i] = i_rc
    no_minus = list(contigs_d.values())
    # return contigs, contigs
    return no_minus, minus_contig


def get_degree(up_dict, down_dict):
    degree = {}
    for up in up_dict:
        degree[up] = [0, len(up_dict[up])]
    for down in down_dict:
        if down in degree:
            degree[down] = [len(down_dict[down]), degree[down][1]]
        else:
            degree[down] = [len(down_dict[down]), 0]
    for i in degree:
        dot.node(i, xlabel=str(degree[i]))
    return degree


def clean_link(overlap):
    """
    Remove transitively-inferrible edges.
    Remove orphan contigs ("island").
    Remove non-branching stretches ("tips").
    Assume "bubble" does not exist.
    Arg:
        overlap(list(blast_result)): link info of contigs
    Return:
        cleaned_link(list(blast_result)): clean link
    """
    raw = {(i[0], i[1]): i for i in overlap}
    up_down = defaultdict(set)
    down_up = defaultdict(set)
    for i in overlap:
        up_down[i[0]].add(i[1])
        down_up[i[1]].add(i[0])
    # Remove transitively-inferible edges that across one contig
    # Use while loop to remove all that kinds of edges?
    shortcuts = set()
    # Remove non-branching stretches
    # One step is enough, if longer, may be alternative path
    tips = set()
    for up, down in up_down.items():
        if len(down) == 1:
            continue
        down_down = set()
        for d in down:
            if d in up_down:
                down_down.update(up_down[d])
            else:
                tips.add((up, d))
        shortcuts.update({(up, i) for i in down_down & down})
    # shortcuts between two circles
    between = []
    for down, up in down_up.items():
        n = len(down)
        if len(up) == 1:
            continue
        for u in up:
            # type II
            if u not in down_up:
                tips.add((u, down))
                n -= 1
            else:
                between.append([u, down])
        if n != 1:
            print('multi', up, down)
    between_d = {}
    shortcuts_b = set()
    for i in between:
        i2 = [j.replace('_RC_', '') for j in i]
        key = '{}---{}'.format(*sorted(i2))
        print(i, key)
        if key in between_d:
            shortcuts_b.add(tuple(i))
            shortcuts_b.add(tuple(between_d[key]))
        else:
            between_d[key] = tuple(i)
    print('between_s', shortcuts_b)
    cleaned_link = [raw[i] for i in raw if (i not in shortcuts and i not in
                                            tips and i not in shortcuts_b)]
    for i in shortcuts:
        dot.edge(*i, color='red')
    for i in tips:
        dot.edge(*i, color='green')
    for i in shortcuts_b:
        dot.edge(*i, color='orange')
    print('all, shortcuts, tips, between_s, clean')
    print(len(overlap), len(shortcuts),  len(tips), len(shortcuts_b), len(cleaned_link))
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
    dot = Digraph(engine='dot', node_attr={'shape': 'cds'})
    overlap = get_overlap(contigs)
    contigs_no_minus, minus_contigs = remove_minus(overlap, contigs)
    overlap_no_minus = get_overlap(contigs_no_minus)
    up_dict = defaultdict(set)
    # down_id: hit
    down_dict = defaultdict(set)
    for i in overlap_no_minus:
        up_dict[i[0]].add(i[1])
        down_dict[i[1]].add(i[0])
    for i in overlap_no_minus:
        dot.node(i[0])
        dot.node(i[1])
        if i[2] == 'plus':
            dot.edge(i[0], i[1], color='#999999')
        else:
            dot.edge(i[0], i[1], color='#999999', style='dashed', dir='both')
    # remove orphan minus
    overlap_no_minus = [i for i in overlap_no_minus if i[2] != 'minus']
    overlap_clean = clean_link(overlap_no_minus)
    links = []
    scaffold = []
    # assume each seq only occurs once
    print()
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
        circle = (link[0][0] == link[-1][1])
        # qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen, qstart,
        # qend, sstart, send
        seq = contigs_d[link[0][0]]
        for node in link[:-1]:
            down = contigs_d[node[1]]
            seq += down[node[11]:]
        # tail do not have link after itself
        if not circle:
            tail_seq = contigs_d[link[-1][1]]
            seq += tail_seq[link[-1][11]:]
        else:
            seq = seq[:-link[-1][11]]
        seq.id = f'Merged_sequence {len(seq)}bp'
        print(circle, 134502, seq.id)
        # seems circle is ok, non-circle is always bad result
        if circle:
            merged.append(seq)
        else:
            continue
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