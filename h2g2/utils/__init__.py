# pip import
import vcf
from Bio import SeqIO

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
    not_block_variant = list()
    for variant in sorted(variants, key=lambda x: x[0].start):
        if len(block) == 0:
            block.append(variant)
            continue
        
        if block[-1][0].CHROM == variant[0].CHROM and block[-1][0].end > variant[0].start:
            block.append(variant)
        else:
            if len(block) > 1:
                list_block.append(block)
            else:
                not_block_variant.append(block[0])
            block = list()
            block.append(variant)

    if len(block) > 1:
        list_block.append(block)

    return not_block_variant, list_block
