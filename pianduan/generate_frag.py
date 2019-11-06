#!/usr/bin/python3

from Bio import SeqIO


# oryza sativa cp's IR 20k
def s(overlap=500, step=40000):
    result = []
    i = 0
    while i < len(a):
        seq = a[i:i+step]
        seq.id = f'{i}-{i+step}'
        result.append(seq)
        i = i + (step-overlap)
    return result


a = SeqIO.read('./Option_1_SRR1328237-rbcL.fasta', 'fasta')
a = a+a[:20000]
Short = s(overlap=100)
SeqIO.write(Short, 'Short.fasta', 'fasta')
Long = s(overlap=2000)
SeqIO.write(Long, 'Long.fasta', 'fasta')
