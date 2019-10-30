basic_path = "/scratch/bioinf/projects/dga-testing-2/output/diploid_assembly/strandseq/freebayes_GQ100_DP50/HG00733_sra_pbsq1-clr_clustV2-100kb/HG00733_hgsvc_pbsq2-ccs_1000/HG00733_1kg_il25k-npe_sseq/consensus"

rule mapping:
    input:
        query  = basic_path+"/{filename}.h1-un.fasta",
        target = basic_path+"/{filename}.h2-un.fasta"
    output:
        "mapping/{filename}_h1_h2.paf"
    shell:
        "minimap2 -cx asm5 --cs {input.query} {input.target} | sort -k6,6 -k8,8n > {output}"
    
rule vcf:
    input:
        "mapping/{filename}_h1_h2.paf"
    output:
        "variant/{filename}.vcf"
    shell:
        "paftools.js call {input} > {output}"


rule all:
    input:
        "variant/HG00733_sra_pbsq1-clr_1000.vcf"
