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


def add_rc(contigs):
    """
    Add reverse-complement contig
    Args:
        contigs(list(SeqRecord)): contigs
    Return:
        raw_and_rc(list(SeqRecord)): minus contigs
    """
    raw_and_rc = []
    for i in contigs:
        i_rc = i.reverse_complement(id=PREFIX+i.id)
        raw_and_rc.extend([i, i_rc])
    return raw_and_rc


def reverse_link(link):
    new = []
    for i in reversed(link):
        if i.startswith(PREFIX):
            i = i.replace(PREFIX, '')
        else:
            i = PREFIX + i
        new.append(i)
    return tuple(new)


def clean_overlap(overlap):
    """
    Remove transitively-inferrible edges ("shortcuts").
    Remove edges across two circle ("shortcuts_b").
    Remove orphan contigs ("island").
    Remove non-branching stretches ("tips").
    Remove alternative path ("bubble").
    Arg:
        overlap(list(blast_result)): link info of contigs
    Return:
        cleaned_overlap(list(blast_result)): clean link
    """
    def get_dict(overlap_list):
        up_down = defaultdict(set)
        down_up = defaultdict(set)
        for i in overlap_list:
            up_down[i[0]].add(i[1])
            down_up[i[1]].add(i[0])
        return up_down, down_up

    raw = {(i[0], i[1]): i for i in overlap}
    up_down, down_up = get_dict(overlap)
    # Remove transitively-inferible edges that across one contig
    # Use while loop to remove all that kinds of edges?
    shortcuts = set()
    # Remove non-branching stretches
    # One step is enough, if longer, may be alternative path
    tips_u_d = set()
    tips_d_u = set()
    for up, down in up_down.items():
        if len(down) == 1:
            continue
        down_down = set()
        for d in down:
            if d in up_down:
                # if len(down_up[d]) > 1
                down_down.update(up_down[d])
            else:
                tips_u_d.add((up, d))
        shortcuts.update({(up, i) for i in down_down & down})
    # shortcuts between two circles
    between = set()
    for down, up in down_up.items():
        if len(up) == 1:
            continue
        for u in up:
            between.add((u, down))
    shortcuts_b = set()
    for i in between:
        if reverse_link(i) in between:
            shortcuts_b.add(i)
            shortcuts_b.add(reverse_link(i))
    # another direction
    between2 = set()
    for up, down in up_down.items():
        if len(down) == 1:
            continue
        for d in down:
            between.add((up, d))
    for i in between2:
        if reverse_link(i) in between2:
            shortcuts_b.add(i)
            shortcuts_b.add(reverse_link(i))
    # non-branching stretches
    # One step is enough, if longer, may be alternative path
    #tips_u_d = set()
    #tips_d_u = set()
    #for up, down in up_down.items():
    #    if len(down) == 1:
    #        continue
    #    for d in down:
    #        if d not in up_down:
    #            tips_u_d.add((up, d))
    #for down, up in down_up.items():
    #    if len(up) == 1:
    #        continue
    #    for u in up:
    #        if u not in down_up:
    #            tips_d_u.add((u, down))
    #short_tips = tips_u_d | tips_d_u
    #short_tips_r = {reverse_link(i) for i in short_tips}
    #short_tips = short_tips.union(short_tips_r)
    # transitively-inferible edges that across one contig
    # Use while loop to remove all that kinds of edges?
    shortcuts = set()
    for up, down in up_down.items():
        if len(down) == 1:
            continue
        down_down = set()
        for d in down:
            if d in up_down:
                down_down.update(up_down[d])
        shortcuts.update({(up, i) for i in down_down & down})
    exclude = shortcuts & shortcuts_b
    shortcuts = shortcuts - exclude
    shortcuts_b = shortcuts_b - exclude
    # to_remove = shortcuts | short_tips | shortcuts_b
    to_remove = shortcuts |  shortcuts_b
    cleaned_overlap = [raw[i] for i in raw if i not in to_remove]
    # a-b-c, a-d-c
    # or a-b-c-d, a-e
    bubble_path = []
    # long tip, a-b-c-d-a, b-d
    tips = set()
    # other types may ganrao bubble detection
    up_down, down_up = get_dict(cleaned_overlap)
    depth = len(up_down) - 1
    for up, down in up_down.items():
        if len(down) <= 1:
            continue
        path = []
        for d in down:
            step = 0
            p = [up, d]
            while step <= depth:
                # more than one upstream
                if len(down_up[d]) > 1:
                    break
                # tail of linear
                if d not in up_down:
                    break
                if up == d:
                    good = tuple(p[:2])
                    bad = {(up, i) for i in up_down[d]}
                    bad.remove(good)
                    r = {reverse_link(i) for i in bad}
                    tips.update(bad)
                    tips.update(r)
                    break
                if len(up_down[d]) > 1:
                    if up in up_down[d]:
                        t = set(up_down[d])
                        t.remove(up)
                        print('t', t)
                        if t is not None and len(t) != 0:
                            tips.update(set((d, dd) for dd in t))
                            r = set((reverse_link((d, dd)) for dd in t))
                            tips.update(r)
                            print('r', r)
                            print('found tip', tips)
                    break
                d = set(up_down[d])
                d = d.pop()
                p.append(d)
                step += 1
            path.append(p)
        n_tail = len(set([p[-1] for p in path]))
        n_length = len(set([len(p) for p in path]))
        if len(path) == 0:
            continue
        if n_tail == 1:
            if n_length == 1:
                # same length, same head/tail
                print('bubble', path)
                # if bubble, keep the first
                bubble_path.extend(path[1:])
            else:
                print('shortcut', path)
    # tips in another direction
    for down, up in down_up.items():
        if len(up) <= 1:
            continue
        for u in up:
            step = 0
            p = [down, u]
            while step <= depth:
                # more than one downstream
                if len(up_down[u]) > 1:
                    break
                # head of linear
                if u not in down_up:
                    break
                if down == u:
                    good = (p[1], p[0])
                    bad = {(i, down) for i in down_up[u]}
                    bad.remove(good)
                    r = {reverse_link(i) for i in bad}
                    tips.update(bad)
                    tips.update(r)
                u = set(down_up[u])
                u = u.pop()
                p.append(u)
                step += 1
    bubble = set()
    for path in bubble_path:
        for i in range(len(path)-1):
            bubble.add((path[i], path[i+1]))
    print('tips', tips)
    to_remove = to_remove.union(bubble).union(tips)
    to_remove = to_remove.union(bubble)
    cleaned_overlap = [raw[i] for i in raw if i not in to_remove]
    for i in shortcuts:
        dot.edge(*i, color='red')
    for i in tips:
        dot.edge(*i, color='green')
    #for i in short_tips:
    #    continue
    #    dot.edge(*i, color='#88ff88')
    for i in shortcuts_b:
        dot.edge(*i, color='orange')
    for i in bubble:
        dot.edge(*i, color='purple')
    with dot.subgraph(name='legend', node_attr={'shape': 'box'}) as s:
        s.node('shortcuts', color='red', style='filled')
        s.node('tips', color='green', style='filled')
        s.node('shortcuts_b', color='orange', style='filled')
        s.node('bubble', color='purple', style='filled')
        s.node('short_tips', color='#88ff88', style='filled')
    print('all, shortcuts, tips, between_s, bubble, clean')
    print(len(overlap), len(shortcuts),  len(tips), len(shortcuts_b),
          len(bubble), len(cleaned_overlap))
    #print('exclude ', exclude)
    #print('shortcuts', shortcuts)
    #print('tips', tips)
    #print('shortcuts_b', shortcuts_b)
    #print('bubble', bubble)
    return cleaned_overlap


def get_link(contigs):
    """
    Get overlap between contigs.
    Remove overlap that may cause chaos.
    Args:
        contigs(list(SeqRecord)): contigs
    Return:
        contigs_and_rc(list(SeqRecord)): contigs with their reverse-complements
        link(list(blast_result)): link info of contigs
    """
    global dot
    dot = Digraph(engine='dot', node_attr={'shape': 'cds'})
    contigs_and_rc = add_rc(contigs)
    overlap = get_overlap(contigs_and_rc)
    for i in overlap:
        dot.node(i[0])
        dot.node(i[1])
        if i[2] == 'plus':
            dot.edge(i[0], i[1], color='#999999')
        else:
            continue
            print()
            dot.edge(i[0], i[1], color='#999999', style='dashed', dir='both')
    # remove minus
    overlap_no_minus = [i for i in overlap if i[2] != 'minus']
    overlap_clean = clean_overlap(overlap_no_minus)
    # in case of missing bubble
    # twice is enough
    overlap_clean = clean_overlap(overlap_clean)
    up_dict = defaultdict(set)
    down_dict = defaultdict(set)
    for i in overlap_clean:
        up_dict[i[0]].add(i[1])
        down_dict[i[1]].add(i[0])
    print('count')
    print('up-down', [i for i in up_dict.items() if len(i[1]) !=1])
    print('down-up', [i for i in down_dict.items() if len(i[1]) !=1])
    print('count')
    links = []
    scaffold = []
    overlap_dict = {i[0]: i for i in overlap_clean}
    up_dict = {i[0]: i for i in overlap_clean}
    # up/down dict should be same
    down_dict = {i[1]: i for i in up_dict.values()}
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
    for i in links:
        print('->'.join([str((j[0], j[1])) for j in i]))
    return contigs_and_rc, links


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
    global PREFIX
    PREFIX = '_RC_'
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
            record.id = f'{fasta.stem.replace(".fasta", "")}-{idx}'
            record.description = ''
            contigs.append(record)
    shuffle(contigs)
    contigs_no_minus, links = get_link(contigs)
    merged = merge_seq(contigs_no_minus, links)
    SeqIO.write(merged, Path(argv[1]).with_suffix('.merge'), 'fasta')


if __name__ == '__main__':
    merge_contigs()
