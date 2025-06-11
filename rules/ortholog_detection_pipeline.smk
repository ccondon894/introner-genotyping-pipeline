import os

GENOME_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies/"
OUTPUT_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph_blast_results"
GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"
FIGURES = "/scratch1/chris/introner-genotyping-pipeline/figures"
SAMPLES = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749", 
           "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", 
           "RCC692", "RCC693", "RCC833"]



def get_targets(query):
    return [s for s in SAMPLES if s != query]

rule all:
    input:
        # Generate all possible query-target combinations excluding self-comparisons
        [os.path.join(OUTPUT_DIR, f"{query}_{target}.left_flank.sorted.bam")
         for query in SAMPLES
         for target in SAMPLES if target != query] +
        [os.path.join(OUTPUT_DIR, f"{query}_{target}.right_flank.sorted.bam")
         for query in SAMPLES
         for target in SAMPLES if target != query],
        os.path.join(GT_DIR, "genotype_matrix.tsv"),
        os.path.join(GT_DIR, "genotype_matrix_oriented.tsv"),

       
rule extract_flanks:
    input:
        fa = os.path.join(OUTPUT_DIR, "{sample}.candidate_loci_plus_flanks.filtered.similarity_checked.fa")
    output:
        left = os.path.join(OUTPUT_DIR, "{sample}.left_flanks.fa"),
        right = os.path.join(OUTPUT_DIR, "{sample}.right_flanks.fa")
    shell:
        """
        python scripts/extract_flanks.py {input.fa} {output.left} {output.right}
        """

rule bwa_index:
    input:
        genome = os.path.join(GENOME_DIR, "{sample}.vg_paths.fa")
    output:
        # BWA creates multiple index files - we'll use .bwt as a representative
        index = os.path.join(GENOME_DIR, "{sample}.vg_paths.fa.bwt")
    log:
        os.path.join(GENOME_DIR, "logs", "{sample}.bwa_index.log")
    threads: 1  # BWA index is single-threaded
    shell:
        """
        bwa index {input.genome} 2> {log}
        """

rule bwa_align_flanks:
    input:
        query_left = os.path.join(OUTPUT_DIR, "{query}.left_flanks.fa"),
        query_right = os.path.join(OUTPUT_DIR, "{query}.right_flanks.fa"),
        target_genome = os.path.join(GENOME_DIR, "{target}.vg_paths.fa"),
        target_index = os.path.join(GENOME_DIR, "{target}.vg_paths.fa.bwt")
    output:
        left_bam = os.path.join(OUTPUT_DIR, "{query}_{target}.left_flank.sorted.bam"),
        right_bam = os.path.join(OUTPUT_DIR, "{query}_{target}.right_flank.sorted.bam"),
        left_bai = os.path.join(OUTPUT_DIR, "{query}_{target}.left_flank.sorted.bam.bai"),
        right_bai = os.path.join(OUTPUT_DIR, "{query}_{target}.right_flank.sorted.bam.bai")
    threads: 4
    shell:
        """
        # Align left flanks
        bwa mem -t {threads} {input.target_genome} {input.query_left} | \
        samtools sort -@ {threads} -o {output.left_bam}
        samtools index {output.left_bam}

        # Align right flanks
        bwa mem -t {threads} {input.target_genome} {input.query_right} | \
        samtools sort -@ {threads} -o {output.right_bam}
        samtools index {output.right_bam}
        """

rule pair_orthologs:
    input:
        left_bam = [os.path.join(OUTPUT_DIR, f"{query}_{target}.left_flank.sorted.bam")
         for query in SAMPLES
         for target in SAMPLES if target != query],
        right_bam = [os.path.join(OUTPUT_DIR, f"{query}_{target}.right_flank.sorted.bam")
         for query in SAMPLES
         for target in SAMPLES if target != query],
    output:
        orthologs = os.path.join(GT_DIR, "ortholog_results.tsv")
    shell:
        """
        python scripts/introner_genotype_matrix_builder.py {OUTPUT_DIR} {output.orthologs}
        """

rule build_genotype_matrix:
    input:
        orthologs = os.path.join(GT_DIR, "ortholog_results.tsv"),
    output:
        genotype_matrix = os.path.join(GT_DIR, "genotype_matrix.tsv"),
    shell:
        """
        python scripts/introner_network.py {input.orthologs} {output.genotype_matrix}
        """

rule fix_orientations:
    input:
        genotype_matrix = os.path.join(GT_DIR, "genotype_matrix.tsv"),
    output:
        fixed_matrix = os.path.join(GT_DIR, "genotype_matrix_oriented.tsv"),
    params:
        genome_dir = GENOME_DIR
    shell:
        """
        python scripts/check_orientation.py {input.genotype_matrix} {params.genome_dir} {output.fixed_matrix} --threads 20
        """