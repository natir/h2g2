#!/usr/bin/env python

# standard import
import sys
from collections import defaultdict

# pip import
import pysam
import vcf
from Bio import SeqIO

def split(args):

    output = open(args.output, 'w')
    tig2seq = read_ref(args.references)

    variants, vcf2haploid, haploid2vcf = read_vcfs(args.vcf_files)
    
    blocks = generate_block(variants)
    
    print("##fileformat=VCFv4.2", file=output)
    print("##FORMAT=<ID=GT,Number={},Type=String,Description=\"Genotype\">".format(len(vcf2haploid)), file=output)
    
    print("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", *args.vcf_files]), file=output)

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
        hetero = [""] * len(args.vcf_files)
        print(alts)
        print()
        for i, alt in enumerate(alts):
            print(i, alt)
            if alt != "":
                if alt not in alts2id:
                    alts2id[alt] = i + 1
                    altid = alts2id[alt]
                    variant = block[i][0]
                    print("variant", variant)
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

def dist(a, b):
    if a > b:
        return a - b
    else:
        return b - a 

        
def read_ref(path):
    tig2seq = dict()
    with open(path, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            tig2seq[record.id] = record.seq

    return tig2seq

def read_vcfs(paths):
    variants = list()
    vcf2haploid = dict()
    haploid2vcf = dict()
    haplo_id = 0
    for vcf_file in paths:
        haplo_id += 1
        vcf_reader = vcf.Reader(open(vcf_file))

        vcf2haploid[vcf_file] = haplo_id
        haploid2vcf[haplo_id] = vcf_file
        for record in vcf_reader:
            variants.append((record, vcf_file))

    return variants, vcf2haploid, haploid2vcf

def generate_block(variants):
    block = list()
    list_block = list()
    for variant in sorted(variants, key=lambda x: x[0].start):
        if len(block) == 0:
            block.append(variant)
            continue
        
        if block[-1][0].CHROM == variant[0].CHROM and block[-1][0].end > variant[0].start:
            block.append(variant)
        else:
            if len(block) > 1:
                list_block.append(block)
            block = list()
            block.append(variant)

    if len(block) > 1:
        list_block.append(block)

    return list_block
