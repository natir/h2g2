#!/usr/bin/env python

# standard import
import sys
from collections import defaultdict

# project import
from ..utils import *

def extract(args):

    output_prefix = args.output
    
    tig2seq = read_ref(args.references)

    variants, vcf2haploid, haploid2vcf = read_vcfs(args.vcf_files)
    
    _, blocks = generate_block(variants)
        
    for block in filter_snp_block(blocks):
        poss = list()
        for v in block:
            poss.append(v[0].start)
            poss.append(v[0].end)
        poss.sort()
        
        begin = poss[0] - args.kmer_size + 1 
        end = poss[-1] + args.kmer_size + 1
        chrom = block[0][0].CHROM
        ref = str(tig2seq[chrom][begin-1:end-1])
    
        alts = {args.references: ref}
        hetero = [""] * len(args.vcf_files)
        for v in block:
            variant = v[0]
            
            change_pos_begin = dist(begin, (variant.start + 1))
            change_pos_end = dist(len(ref), (end - (variant.end + 1)))
            alt = str(variant.ALT[0])
            alt = ref[:change_pos_begin]+alt+ref[change_pos_end:]
            alts[v[1]] = alt

        output = open(output_prefix+chrom+'_'+str(begin)+'.fasta', 'w')
        print("\n".join([">{}\n{}".format(id, seq) for id, seq in alts.items()]), file=output)


def filter_snp_block(blocks):
    for block in blocks:
        if all(v[0].is_snp for v in block):
            yield block
            
