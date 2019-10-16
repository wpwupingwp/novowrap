#!/usr/bin/python3

from Bio import SeqIO


def s(overlap=500, step=50000):
    result = []
    i = 0
    while i < len(a):
        seq = a[i:i+step]
        seq.id = f'{i}-{i+step}'
        result.append(seq)
        i = i + (step-overlap)
    return result


a = SeqIO.read('./Option_1_SRR1328237-rbcL.fasta', 'fasta')
a500b = s(overlap=500, step=20000)
SeqIO.write(a500b, 'Oryza_sativa_500bp.fasta', 'fasta')
a1k = s(overlap=1000, step=20000)
SeqIO.write(a1k, 'Oryza_sativa_1k.fasta', 'fasta')
a2k = s(overlap=2000, step=20000)
SeqIO.write(a2k, 'Oryza_sativa_2k.fasta', 'fasta')
a5k = s(overlap=5000, step=20000)
SeqIO.write(a5k, 'Oryza_sativa_5k.fasta', 'fasta')
