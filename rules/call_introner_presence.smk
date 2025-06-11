import os

GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"
OUTPUT_DIR = "/scratch1/chris/introner-genotyping-pipeline/bam_coverage_output"
BAM_DIR = "/storage1/chris/introner-genotyping-bams"
GENOME_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies/"
READS_DIR = "/scratch1/alex/pusilla/data/reads"
FIGURES = "/scratch1/chris/introner-genotyping-pipeline/figures"


strains = ["RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
paths = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]

rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{path}.loci.bed"), path=paths),
        expand(os.path.join(OUTPUT_DIR, "{strain}.locus_depth.tsv"), strain=strains),
        expand(os.path.join(OUTPUT_DIR, "{strain}.loci.coverage_calls.bed"), strain=strains),
        os.path.join(GT_DIR, "genotype_matrix_updated.tsv"),
        os.path.join(GT_DIR, "genotype_matrix_summary.txt"),

rule align_reads:
    input:
        fa = os.path.join(GENOME_DIR, "{strain}.vg_paths.fa"),
        fq1 = os.path.join(READS_DIR, "{strain}/{strain}.R1.sorted.fastq"),
        fq2 = os.path.join(READS_DIR, "{strain}/{strain}.R2.sorted.fastq")
    output:
        bam = os.path.join(BAM_DIR, "{strain}.sorted.bam")
    threads: 8
    shell:
        """
        bwa index {input.fa} 2> /dev/null
        bwa mem -t {threads} {input.fa} {input.fq1} {input.fq2} | \
        samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        """

rule make_loci_beds:
    input:
        gt_matrix = os.path.join(GT_DIR, "genotype_matrix_oriented.tsv"),
    
    output:
        bed = os.path.join(OUTPUT_DIR, "{path}.loci.bed")
    
    shell:
        """
        python scripts/make_sample_bed_from_gt_matrix.py {input.gt_matrix} {wildcards.path} {output.bed}
        """

rule samtools_depth:
    input:
        bed = os.path.join(OUTPUT_DIR, "{strain}.loci.bed"),
        bam = os.path.join(BAM_DIR, "{strain}.sorted.bam")
    output:
        os.path.join(OUTPUT_DIR, "{strain}.locus_depth.tsv")
    shell:
        """
        # Skip header line (if present) by using awk to process only lines that don't start with #, track or header
        awk '!/^#/ && !/^track/ && !/^browser/ && !/^contig/' {input.bed} | samtools depth -Q 30 -a -J -b - {input.bam} > {output}
        """

rule call_introner_presence:
    input:
        bed = os.path.join(OUTPUT_DIR, "{strain}.loci.bed"),
        tsv = os.path.join(OUTPUT_DIR, "{strain}.locus_depth.tsv"),
        gt_matrix = os.path.join(GT_DIR, "genotype_matrix_oriented.tsv")
    output:
        bed = os.path.join(OUTPUT_DIR, "{strain}.loci.coverage_calls.bed"),

    log: "logs/{strain}.coverage_calls.log"
    shell:
        """
        python scripts/call_locus_from_bam_coverage_v2.py {input.bed} {input.tsv} {input.gt_matrix} {output.bed} > {log}
        """

rule update_gt_matrix:
    input:
        gt_matrix = os.path.join(GT_DIR, "genotype_matrix_oriented.tsv"),
        beds = expand(os.path.join(OUTPUT_DIR, "{strain}.loci.coverage_calls.bed"), strain=strains)
    output:
       gt_matrix = os.path.join(GT_DIR, "genotype_matrix_updated.tsv"),
    params:
        bed_path = OUTPUT_DIR
    shell:
        """
        python scripts/update_gt_matrixes_bam_cov.py {input.gt_matrix} {params.bed_path} {output.gt_matrix}
        """

rule plot_genotype_matrix_data:
    input:
        genotype_matrix = os.path.join(GT_DIR, "genotype_matrix_updated.tsv"),
    output:
        plot1 = os.path.join(FIGURES, f"genotype_matrix_analysis_28052025_2d_afs.png"),
        summary = os.path.join(GT_DIR, "genotype_matrix_summary.txt"),
    params:
        prefix = "genotype_matrix_analysis_28052025"
    shell:
        """
        python scripts/genotype_matrix_data_viz.py {input.genotype_matrix} {params.prefix} > {output.summary}
        """
