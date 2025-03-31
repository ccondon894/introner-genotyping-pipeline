import os

GENOME_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies/"
OUTPUT_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph_blast_results"
GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"
SAMPLES = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", 
           "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", 
           "RCC692", "RCC693", "RCC833", "RCC835"]

LOCUS_THRESHOLD = 0.8  # Minimum sequence similarity for locus presence

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
        # [os.path.join(OUTPUT_DIR, f"{query}.ortholog_results.adjusted.tsv")
        #  for query in SAMPLES] +
        # Add index files
        os.path.join(GT_DIR, "genotype_matrix.metadata.4.tsv"),
        os.path.join(GT_DIR, "genotype_matrix.calls.4.tsv"),

rule extract_flanks:
    input:
        fa = os.path.join(OUTPUT_DIR, "{sample}.candidate_loci_plus_flanks.filtered.fa")
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

rule analyze_bam_alignments:
    input:
        left_bam = os.path.join(OUTPUT_DIR, "{query}_{target}.left_flank.sorted.bam"),
        right_bam = os.path.join(OUTPUT_DIR, "{query}_{target}.right_flank.sorted.bam"),
        query_fa = os.path.join(OUTPUT_DIR, "{query}.candidate_loci_plus_flanks.filtered.fa"),
        target_fa = os.path.join(GENOME_DIR, "{target}.vg_paths.fa")
    output:
        tsv = os.path.join(OUTPUT_DIR, "{query}_{target}.ortholog_pairs.tsv")
    shell:
        """
        python scripts/analyze_bam_alignments.py \
            --left_bam {input.left_bam} \
            --right_bam {input.right_bam} \
            --query_fa {input.query_fa} \
            --target_fa {input.target_fa} \
            --output {output.tsv}
        """
        

rule combine_results:
    input:
        ortholog_pairs = lambda wildcards: expand(os.path.join(OUTPUT_DIR, "{{query}}_{target}.ortholog_pairs.tsv"),
                              target=[s for s in SAMPLES if s != wildcards.query])
    output:
        combined = os.path.join(OUTPUT_DIR, "{query}.ortholog_results.tsv")
    shell:
        """
        python scripts/combine_ortholog_results.py \
            --input_files {input.ortholog_pairs} \
            --output {output.combined}
        """

rule adjust_coordinates:
    input:
        combined = os.path.join(OUTPUT_DIR, "{query}.ortholog_results.tsv"),
        bed_dir = OUTPUT_DIR
    output:
        combined = os.path.join(OUTPUT_DIR, "{query}.ortholog_results.adjusted.tsv"),
    shell:
        """
        python scripts/adjust_bed_coordinates.py {input.combined} {input.bed_dir} {output.combined}
        """

rule combine_all_results:
    input:
        sample_results = expand(os.path.join(OUTPUT_DIR, "{query}.ortholog_results.adjusted.tsv"), query=SAMPLES),
    output:
        combined = os.path.join(OUTPUT_DIR, "combined_orthologs.tsv"),
    shell:
        """
        python scripts/make_initial_genotype_matrix.py {output.combined}
        """

rule make_genotype_matrixes:
    input:
        combined = os.path.join(OUTPUT_DIR, "combined_orthologs.tsv"),
    output:
        md = os.path.join(GT_DIR, "genotype_matrix.metadata.4.tsv"),
        calls = os.path.join(GT_DIR, "genotype_matrix.calls.4.tsv")
    shell:
        """
        python scripts/make_genotype_matrix_v2.py {input.combined} {output.md} {output.calls}
        """