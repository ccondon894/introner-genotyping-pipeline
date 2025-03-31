import os

GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"
OUTPUT_DIR = "/scratch1/chris/introner-genotyping-pipeline/bam_coverage_output"
BAM_DIR = "/storage1/chris/introner-genotyping-bams"
GENOME_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies/"
READS_DIR = "/scratch1/alex/pusilla/data/reads"

strains = ["CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]
paths = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{path}.loci.bed"), path=paths),
        expand(os.path.join(OUTPUT_DIR, "{strain}.locus_depth.tsv"), strain=strains),
        expand(os.path.join(OUTPUT_DIR, "{strain}.loci.coverage_calls.bed"), strain=strains),
        os.path.join(GT_DIR, "genotype_matrix.metadata.4.updated.tsv"),
        os.path.join(GT_DIR, "genotype_matrix.calls.4.updated.tsv"),

rule align_reads:
    input:
        fa = os.path.join(GENOME_DIR, "{strain}.vg_paths.fa"),
        fq1 = os.path.join(READS_DIR, "{strain}/{strain}.R1.sorted.fastq"),
        fq2 = os.path.join(READS_DIR, "{strain}/{strain}.R2.sorted.fastq")
    output:
        bam = os.path.join(BAM_DIR, "{strain}.sorted.bam")
    shell:
        """
        bwa index {input.fa} 2> /dev/null
        bwa mem -t 8 {input.fa} {input.fq1} {input.fq2} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule make_loci_beds:
    input:
        metadata = os.path.join(GT_DIR, "genotype_matrix.metadata.4.tsv"),
    
    output:
        bed = os.path.join(OUTPUT_DIR, "{path}.loci.bed")
    
    shell:
        """
        python scripts/make_loci_beds_from_metadata.py {input.metadata} {wildcards.path} {output.bed}
        """

rule samtools_depth:
    input:
        bed = os.path.join(OUTPUT_DIR, "{strain}.loci.bed"),
        bam = os.path.join(BAM_DIR, "{strain}.sorted.bam")
    output:
        os.path.join(OUTPUT_DIR, "{strain}.locus_depth.tsv")
    shell:
        """
        samtools depth -Q 30 -a -J -b {input.bed} {input.bam} > {output}
        """

rule call_introner_presence:
    input:
        bed = os.path.join(OUTPUT_DIR, "{strain}.loci.bed"),
        tsv = os.path.join(OUTPUT_DIR, "{strain}.locus_depth.tsv"),
        call_file = os.path.join(GT_DIR, "genotype_matrix.calls.4.tsv")
    output:
        bed = os.path.join(OUTPUT_DIR, "{strain}.loci.coverage_calls.bed"),

    log: "logs/{strain}.coverage_calls.log"
    shell:
        """
        python scripts/call_locus_from_bam_coverage.py {input.bed} {input.tsv} {input.call_file} {output.bed} > {log}
        """

rule update_gt_matrixes:
    input:
        meta = os.path.join(GT_DIR, "genotype_matrix.metadata.4.tsv"),
        call = os.path.join(GT_DIR, "genotype_matrix.calls.4.tsv"),
        beds = expand(os.path.join(OUTPUT_DIR, "{strain}.loci.coverage_calls.bed"), strain=strains)
    output:
        meta = os.path.join(GT_DIR, "genotype_matrix.metadata.4.updated.tsv"),
        call = os.path.join(GT_DIR, "genotype_matrix.calls.4.updated.tsv")
    params:
        bed_path = OUTPUT_DIR
    shell:
        """
        python scripts/update_gt_matrixes_bam_cov.py {input.meta} {input.call} {params.bed_path} {output.meta} {output.call}
        """


