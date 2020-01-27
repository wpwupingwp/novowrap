#!/usr/bin/python3

from collections import defaultdict, deque
from itertools import product as cartesian_product
from pathlib import Path
from random import shuffle
import argparse
import logging

from Bio import SeqIO
try:
    from graphviz import Digraph
    have_dot = True
except ImportError:
    have_dot = False

from utils import blast, parse_blast_tab

PREFIX = '_RC_'

# define logger
FMT = '%(asctime)s %(levelname)-8s %(message)s'
DATEFMT = '%H:%M:%S'
logging.basicConfig(format=FMT, datefmt=DATEFMT, level=logging.INFO)
try:
    import coloredlogs
    coloredlogs.install(level=logging.INFO, fmt=FMT, datefmt=DATEFMT)
except ImportError:
    pass
# inherit logger from novowrap
log = logging.getLogger('novowrap')


def parse_args(arg_list=None):
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('input', nargs='+', help='input filenames')
    arg.add_argument('-o', '-out', dest='out', help='output filename')
    if arg_list is None:
        return arg.parse_args()
    else:
        return arg.parse_args(arg_list)


def reverse_link(link):
    """
    Reverse blast_result
    """
    new = []
    for i in reversed(link):
        if i.startswith(PREFIX):
            i = i.replace(PREFIX, '')
        else:
            i = PREFIX + i
        new.append(i)
    return tuple(new)


def get_path(overlap):
    """
    Get path from overlap information.
    Args:
        overlap(list(blast_result)): overlaps
    Yield:
        path(list(blast_result)): ordered path
        is_circle(bool): is circle or not (linear)
    """
    scaffold = deque()
    overlap_dict = {i[0]: i for i in overlap}
    up_dict = {i[0]: i for i in overlap}
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
            up_name = None
            down_name = None
        # if downstream not in overlap_dict:
        if up_name is None:
            pass
        elif up_name in down_dict:
            upstream = down_dict[up_name][0]
            if upstream in overlap_dict:
                scaffold.appendleft(overlap_dict.pop(upstream))
            else:
                scaffold.appendleft([None, None])
        else:
            # [None, None] as head/tail
            scaffold.appendleft([None, None])
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
            scaffold = list(scaffold)[1:-1]
            is_circle = (scaffold[0][0] == scaffold[-1][1])
            yield scaffold, is_circle
            scaffold = deque()
            try:
                scaffold.append(overlap_dict.popitem()[1])
            except KeyError:
                break


def get_contig(files, out):
    """
    Get contigs (SeqRecord) from input files.
    Add reverse-complement of each contigs.
    Args:
        files(list(str)): input files
        out(Path): output file name
    Return:
        contigs_and_rc(list(SeqRecord)): contigs with their reverse-complements
        contigs_and_rc_fasta(Path): fasta file
    """
    contigs_and_rc = []
    for f in files:
        fasta = Path(f)
        # some id may already used PREFIX
        fasta_stem = fasta.stem.replace('.fasta', '').replace(PREFIX, '-RC-')
        for idx, record in enumerate(SeqIO.parse(fasta, 'fasta')):
            record.id = f'{fasta_stem}-{idx}'
            record.description = ''
            contigs_and_rc.append(record)
            record_rc = record.reverse_complement(id=PREFIX+record.id)
            contigs_and_rc.append(record_rc)
    # print(necessary?)
    shuffle(contigs_and_rc)
    log.info(f'Got {len(contigs_and_rc)//2} contigs from input files.')
    contigs_and_rc_fasta = out.with_suffix('.with_rc.fasta')
    SeqIO.write(contigs_and_rc, contigs_and_rc_fasta, 'fasta')
    return contigs_and_rc, contigs_and_rc_fasta


def get_overlap(contigs, contigs_and_rc_fasta):
    """
    Get overlap by BLAST.
    Args:
        contigs_and_rc_fasta(Path): fasta file
        contigs(list(SeqRecord)): contigs
    Return:
        link_info(list(blast_result)): link info of contigs
    """
    blast_result, blast_log = blast(contigs_and_rc_fasta, contigs_and_rc_fasta)
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
            # normally, subject must be object's downstream
            # for tail->head, may be partial match
            # ambiguous base should not at head/tail
            if sstrand == 'plus':
                if sstart == 1 and qend == qlen:
                    pass
                else:
                    continue
            else:
                # skip minus
                continue
            overlap.append(hit)
    return overlap


def clean_overlap(overlap):
    """
    Remove transitively-inferrible edges ("shortcuts").
    Remove edges across two circle ("shortcuts_b").
    Remove non-branching stretches ("tips").
    Remove alternative path ("bubble").
    Arg:
        overlap(list(blast_result)): link info of contigs
    Return:
        cleaned_overlap(list(blast_result)): clean link
        edges(dict(set(dot.edge))): edges of overlap
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
    # shortcuts_b in another direction
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
    # Remove transitively-inferible edges that across one contig
    # One step is enough, if longer, may be alternative path
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
    to_remove = shortcuts | shortcuts_b
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
                        if t is not None and len(t) != 0:
                            tips.update(set((d, dd) for dd in t))
                            r = set((reverse_link((d, dd)) for dd in t))
                            tips.update(r)
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
        # same length, same head/tail
        if n_tail == 1 and n_length == 1:
            # if bubble, keep the first
            bubble_path.extend(path[1:])
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
    to_remove = to_remove.union(bubble).union(tips)
    to_remove = to_remove.union(bubble)
    cleaned_overlap = [raw[i] for i in raw if i not in to_remove]
    edges = {'shortcuts': shortcuts, 'tips': tips,
             'shortcuts_b': shortcuts_b, 'bubble': bubble,
             'exclude': exclude}
    return cleaned_overlap, edges


def get_link(contigs_and_rc, contigs_and_rc_fasta):
    """
    Get overlap between contigs.
    Remove overlap that may cause chaos.
    Args:
        contigs_and_rc(list(SeqRecord)): contigs with their reverse-complements
        contigs_and_rc_fasta(Path): fasta file
    Return:
        link(list(blast_result)): link info of contigs
    """
    MAX_TRY = 2 ** 16
    overlap = get_overlap(contigs_and_rc, contigs_and_rc_fasta)
    overlap_clean, edges = clean_overlap(overlap)
    if len(overlap_clean) == 0:
        log.critical('Bad input. Please check again.')
        return []
    overlap_clean_dict = {(i[0], i[1]): i for i in overlap_clean}
    for i in edges:
        log.debug(f'{i}, {len(edges[i])}')
    up_dict = defaultdict(set)
    down_dict = defaultdict(set)
    for i in overlap_clean:
        up_dict[i[0]].add(i[1])
        down_dict[i[1]].add(i[0])
    links = []
    possible = []
    up_down_clean = []
    for up, down in up_dict.items():
        if len(down) > 1:
            possible.append([overlap_clean_dict[(up, d)] for d in down])
        else:
            up_down_clean.append(overlap_clean_dict[(up, down.pop())])
    n = 0
    if len(possible) > MAX_TRY:
        log.debug(f'Too many possible ({len(possible)}).')
    for i in cartesian_product(*possible):
        if n >= MAX_TRY:
            break
        n += 1
        combine = up_down_clean + list(i)
        for link in get_path(combine):
            links.append(link)
    edges['link'] = set()
    for i in links:
        for j in i[0]:
            edges['link'].add((j[0], j[1]))
    # draw
    dot_out = contigs_and_rc_fasta.with_suffix('.dot')
    if have_dot:
        dot = Digraph(engine='dot', node_attr={'shape': 'cds'})
        for i in overlap:
            dot.node(i[0])
            dot.node(i[1])
        # dot.edge(i[0], i[1], color='#999999')
        for edge in edges['link']:
            dot.edge(*edge, color='blue')
        for edge in edges['shortcuts']:
            dot.edge(*edge, color='red')
        for edge in edges['tips']:
            dot.edge(*edge, color='green')
        for edge in edges['shortcuts_b']:
            dot.edge(*edge, color='orange')
        for edge in edges['bubble']:
            dot.edge(*edge, color='purple')
        dot.render(dot_out)
    else:
        log.debug('Cannot find graphviz, skip drawing figure.')
    # for i in links: print('->'.join([str((j[0], j[1])) for j in i[0]]))
    return links


def merge_seq(contigs, links):
    """
    Use overlap information of contigs to merge them.
    Assume link_info only contains one kind of link route.
    Arg:
        contigs(list(SeqRecord)): contigs
        links(list(blast_result, is_circle)): link info of contigs
    Return:
        merged(list(SeqRecord)): list of merged sequences
    """
    # give 2 circular result at most
    MAX_CIRCLE = 2
    # break if too much linear
    MAX_LINEAR = 20
    contigs_d = {i.id: i for i in contigs}
    circular = []
    linear_d = {}
    for link, is_circle in links:
        if len(circular) >= MAX_CIRCLE:
            break
        elif len(linear_d) >= MAX_LINEAR:
            log.debug(f'Cannot find circular assembly. Instead found '
                      '{len(circular)} linear asseblies. Break.')
            break
        # qseqid, sseqid, sstrand, qlen, slen, length, pident, gapopen, qstart,
        # qend, sstart, send
        seq = contigs_d[link[0][0]]
        for node in link[:-1]:
            down = contigs_d[node[1]]
            seq += down[node[11]:]
        # tail do not have link after itself
        if not is_circle:
            tail_seq = contigs_d[link[-1][1]]
            try:
                seq += tail_seq[link[-1][11]:]
            # rare
            except IndexError:
                pass
        else:
            seq = seq[:-link[-1][11]]
        seq.id = f'Merged_sequence {len(seq)}bp is_circle={is_circle}'
        seq.description = ''
        # seems circle is ok, non-circle is always bad result
        if is_circle:
            circular.append(seq)
        else:
            linear_d[len(seq)] = seq
    # keep longest
    linear_long = sorted(list(linear_d.values()), key=len,
                         reverse=True)[:MAX_LINEAR]
    if len(circular) != 0:
        return circular
    else:
        log.warning('No circular assembly found.')
        return linear_long


def merge_main(arg_str=None):
    if arg_str is None:
        arg = parse_args()
    else:
        arg = parse_args(arg_str.split(' '))
    if arg.out is None:
        arg.out = Path(arg.input[0]).with_suffix('.merge')
    else:
        arg.out = Path(arg.out)
    contigs_and_rc, contigs_and_rc_fasta = get_contig(arg.input, arg.out)
    links = get_link(contigs_and_rc, contigs_and_rc_fasta)
    if len(links) == 0:
        merged = []
    else:
        merged = merge_seq(contigs_and_rc, links)
    if len(merged) == 0:
        log.critical('Failed to assembly contigs.')
    else:
        log.info(f'Got {len(merged)} assemblies.')
        SeqIO.write(merged, arg.out, 'fasta')
    return len(merged), arg.out


if __name__ == '__main__':
    merge_main()
