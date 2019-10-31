import os

rule mapping:
    input:
        "assembly/asm_{asmid}.fasta"
    output:
        "mapping/ref_asm_{asmid}.sam"
    shell:
        "minimap2 -cx asm5 --cs assembly/reference.fasta {input} | sort -k6,6 -k8,8n > {output}"
    
rule vcf:
    input:
        "mapping/ref_asm_{asmid}.sam"
    output:
        "variant/ref_asm_{asmid}.vcf"
    shell:
        "paftools.js call -f assembly/reference.fasta {input} > {output}"

rule all:
    input:
        ["variant/ref_{}.vcf".format(os.path.splitext(entry.name)[0]) for entry in os.scandir("assembly") if entry.is_file() and entry.name.startswith("asm")]

