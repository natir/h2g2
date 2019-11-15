#!/usr/bin/env python

# standard import
import os
import sys

# pip import
from Bio import SeqIO

def prepare(args):
    if args.reference:
        if len(args.assemblies) < 1:
            print("You need to provide a reference and at least one assembly", file=sys.stderr)
            return 1
        else:
            reference = args.reference
            assemblies = args.assemblies
    else:
        if len(args.assemblies) < 2:
            print("You need to provide more than two assembly", file=sys.stderr)
            return 1
        reference, assemblies = get_ref_by_N50(args.assemblies)
    
    os.symlink(os.path.abspath(reference), os.path.abspath("assembly/reference.fasta"))
    for i, asm in enumerate(assemblies):
        os.symlink(os.path.abspath(asm), os.path.abspath("assembly/asm_{}.fasta".format(i)))

def get_ref_by_N50(assemblies):
    if len(assemblies) < 2:
        print("You need to provide more than two assembly", file=sys.stderr)
        return 1

    filename2N50 = dict()
    
    for asm in assemblies:
        filename2N50[asm] = __get_N(asm, 0.5)
        
    filename_sort = sorted(filename2N50, key=filename2N50.get)

    return filename_sort[0], filename_sort[1:]
        
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
