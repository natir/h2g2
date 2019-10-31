#!/usr/bin/env python

# standard import
import os

# pip import
from Bio import SeqIO

def prepare(args):
    if len(args.assembly) < 2:
        print("You need to provide more than two assembly", file=sys.stderr)
        return 1

    filename2N50 = dict()
    
    for asm in args.assembly:
        filename2N50[asm] = __get_N(asm, 0.5)
        
    filename_sort = sorted(filename2N50, key=filename2N50.get)
        
    os.symlink(os.path.abspath(filename_sort.pop()), os.path.abspath("assembly/reference.fasta"))
    for i, asm in enumerate(filename_sort):
        os.symlink(os.path.abspath(asm), os.path.abspath("assembly/asm_{}.fasta".format(i)))

def __get_N(filename: str, n: float):
    seqlen = list()
    total_len = 0
    
    for record in SeqIO.parse(filename, "fasta"):
        seqlen.append(len(record.seq))
        total_len += len(record.seq)
        
    seqlen.sort()

    threshold = total_len * n
    cumulative_len = 0
    for l in seqlen:
        cumulative_len += l
        if cumulative_len > threshold:
            return l

    return -1
