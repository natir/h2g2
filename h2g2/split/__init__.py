#!/usr/bin/env python

# standard import
import sys
from collections import defaultdict

# project import
from ..utils import *

def split(args):

    output = open(args.output, 'w')
    tig2seq = read_ref(args.references)

    variants, vcf2haploid, haploid2vcf = read_vcfs(args.vcf_files)
    
    others, blocks = generate_block(variants)
    
    print("##fileformat=VCFv4.2", file=output)
    print("##FORMAT=<ID=GT,Number={},Type=String,Description=\"Genotype\">".format(len(vcf2haploid)), file=output)
    
    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *args.vcf_files]), file=output)

    for (v, f) in others:
        gt = ["0|0"] * len(args.vcf_files)
        gt[args.vcf_files.index(f)] = v.genotype(v.samples[0].sample)['GT']

        print("\t".join([v.CHROM, # #CHROM
                         str(v.POS), # POS
                         '.', # ID
                         v.REF, # REF
                         ','.join([s.sequence for s in v.ALT]), # ALT
                         '.', # QUAL
                         '.', # FILTER
                         '.', # INFO
                         'GT', # FORMAT
                         '\t'.join(gt)
                         ]),
              file=output
        )
    
    for block in blocks:
        poss = list()
        for v in block:
            poss.append(v[0].start)
            poss.append(v[0].end)
        poss.sort()
        
        begin = poss[0] + 1
        end = poss[-1] + 1
        chrom = block[0][0].CHROM
        ref = str(tig2seq[chrom][begin-1:end-1])

        alts = [""] * len(args.vcf_files)

        for v in block:
            variant = v[0]

            change_pos_begin = dist(begin, (variant.start + 1))
            change_pos_end = dist(len(ref), (end - (variant.end + 1)))

            alt = str(variant.ALT[0])
            alt = ref[:change_pos_begin]+alt+ref[change_pos_end:]

            alts[vcf2haploid[v[1]] - 1] = alt

        corrected_alts = []
        alts2id = dict()
        hetero = ["0|0"] * len(args.vcf_files)
        i = -1
        for alt in alts:
            if alt != "":
                if alt not in alts2id:
                    i += 1
                    alts2id[alt] = i + 1
                    altid = alts2id[alt]
                    variant = block[i][0]
                    corrected_alts.append(alt)
                else:
                    altid = alts2id[alt]
                hetero[i] = variant.genotype(variant.samples[0].sample)['GT'].replace("1", str(altid))
            else:
                hetero[i] = "0|0"
            
            
        print("\t".join([chrom, # #CHROM
                        str(begin), # POS
                        '.',   # ID
                        ref,   # REF
                        ','.join(corrected_alts), # ALT
                        '.',   # QUAL
                        '.',   # FILTER
                        '.',   # INFO
                        'GT',   # FORMAT
                        '\t'.join(hetero) 
                        ]),
              file=output
        )
