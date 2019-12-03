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
    arg.add_argument('-extend', type=int, default=80000, help='extend length')
    arg.add_argument('-overlap', type=int, default=1000, help='overlap length')
    arg.add_argument('-o', '--out', default='out',
                     help='output directory')
    return arg.parse_args()


def generate(overlap, step):
    result = []
    i = 0
    while i < len(raw):
        seq = raw[i:i+step]
        seq.id = f'{i}-{i+step}'
        result.append(seq)
        i = i + (step-overlap)
    return result


arg = parse_args()
print(arg)
raw = SeqIO.read(arg.file, 'fasta')
raw = raw + raw[:arg.extend]
Short = generate(overlap=arg.overlap//2, step=arg.step)
print('Short:', len(Short))
SeqIO.write(Short, 'Short.fasta', 'fasta')
Long = generate(overlap=arg.overlap, step=arg.step)
print('Long:', len(Long))
SeqIO.write(Long, 'Long.fasta', 'fasta')
