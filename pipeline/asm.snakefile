import os

rule asm:
    input:
        "hard_region/sequences/{chr}_{begin}.fasta"
    output:
        "hard_region/assemblies/{chr}_{begin}_k{k}.unitigs.fa"
    shell:
        "bcalm -in {input} -out hard_region/assemblies/{wildcards.chr}_{wildcards.begin}_k{wildcards.k} -kmer-size {wildcards.k} -abundance-min 0"

rule all:
    input:
        ["hard_region/assemblies/"+os.path.splitext(entry.name)[0]+"_k17.unitigs.fa" for entry in os.scandir("hard_region/sequences/") if entry.is_file()]
