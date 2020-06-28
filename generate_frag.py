#!/usr/bin/python3

from Bio import SeqIO
import argparse


# oryza sativa cp's IR 20k

def parse_args():
    arg = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    arg.add_argument('-file', default='./Option_1_SRR1328237-rbcL.fasta',
                     help='input file')
    arg.add_argument('-step', type=int, default=50000, help='step length')
    arg.add_argument('-extend', type=int, default=40000, help='extend length')
    arg.add_argument('-overlap', type=int, default=1000, help='overlap length')
    return arg.parse_args()


def generate(raw_extend, overlap, step):
    result = []
    i = 0
    while i < len(raw_extend, raw):
        seq = raw[i:i+step]
        seq.id = f'{i}-{i+step}'
        result.append(seq)
        i = i + (step-overlap)
    return result


arg = parse_args()
raw = SeqIO.read(arg.file, 'fasta')
raw_extend = raw + raw[:arg.extend]
contigs = generate(raw_extend, overlap=arg.overlap, step=arg.step)
print('step', arg.step, 'overlap', arg.overlap, 'extend', arg.extend,
      'contigs', len(contigs))
filename = f'L_{arg.step}-O_{arg.overlap}-E_{arg.extend}'
SeqIO.write(contigs, filename, 'fasta')
