import os

# Configuration
GENOME_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies/"
GTF_DIR = "/scratch1/chris/introner-genotyping-pipeline/gene_annotations"
REF_GTF = os.path.join(GTF_DIR, "CCMP1545.vg_paths.gtf")
REF_GENOME = os.path.join(GENOME_DIR, "CCMP1545.vg_paths.fa")
PROTEIN_DIR = "/scratch1/chris/introner-genotyping-pipeline/protein_sequences"  # Directory for extracted protein sequences
STATS_DIR = "/scratch1/chris/introner-genotyping-pipeline/gtf_stats"

# Sample list (excluding reference)
# SAMPLES = ["CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", 
#            "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", 
#            "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]
SAMPLES = ["RCC1749"]

wildcard_constraints:
    sample = "|".join(SAMPLES)

rule all:
    input:
        expand(os.path.join(GTF_DIR, "{sample}.proteins.renamed.gtf"), sample=SAMPLES),
        expand(os.path.join(STATS_DIR, "{sample}_basic_stats.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS_DIR, "{sample}_ref_chrom_dist.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS_DIR, "{sample}_query_chrom_dist.tsv"), sample=SAMPLES),
        expand(os.path.join(STATS_DIR, "{sample}_overlap_stats.tsv"), sample=SAMPLES),
        # os.path.join(STATS_DIR, "combined_stats.tsv")

rule extract_proteins:
    input:
        gtf = REF_GTF,
        genome = REF_GENOME
    output:
        proteins = os.path.join(PROTEIN_DIR, "reference_proteins.fa")
    shell:
        """
        agat_sp_extract_sequences.pl \
            --gff {input.gtf} \
            --fasta {input.genome} \
            --type CDS \
            --merge \
            --protein \
            --output {output.proteins}
        """

rule miniprot_align:
    input:
        proteins = os.path.join(PROTEIN_DIR, "reference_proteins.fa"),
        target_genome = os.path.join(GENOME_DIR, "{sample}.vg_paths.fa")
    output:
        gff = os.path.join(GTF_DIR, "{sample}.raw.gff")
    threads: 8
    shell:
        """
        miniprot -t {threads} --gff {input.target_genome} {input.proteins} > {output.gff}
        """

rule filter_gff:
    input:
        gff = os.path.join(GTF_DIR, "{sample}.raw.gff"),
        ref_gtf = REF_GTF
    output:
        filtered_gff = os.path.join(GTF_DIR, "{sample}.filtered.gff")
    shell:
        """
        python scripts/filter_protein_gtf.py \
            --input {input.gff} \
            --output {output.filtered_gff}
        """

rule convert_to_gtf:
    input:
        gff = os.path.join(GTF_DIR, "{sample}.filtered.gff")
    output:
        gtf = os.path.join(GTF_DIR, "{sample}.proteins.gtf")
    shell:
       """
        agat_convert_sp_gff2gtf.pl \
            --gff {input.gff} \
            --output {output.gtf}
        """

rule rename_gene_ids:
    input:
        gtf = os.path.join(GTF_DIR, "{sample}.proteins.gtf"),
        ref_gtf = REF_GTF
    output:
        gtf = os.path.join(GTF_DIR, "{sample}.proteins.renamed.gtf")
    shell:
        """
        python scripts/rename_gene_ids.py \
        --reference {input.ref_gtf} \
        --query {input.gtf} \
        --output {output.gtf}
        """

rule compute_stats:
    input:
        query_gtf = os.path.join(GTF_DIR, "{sample}.proteins.renamed.gtf"),
        ref_gtf = REF_GTF
    output:
        basic = os.path.join(STATS_DIR, "{sample}_basic_stats.tsv"),
        ref_chrom = os.path.join(STATS_DIR, "{sample}_ref_chrom_dist.tsv"),
        query_chrom = os.path.join(STATS_DIR, "{sample}_query_chrom_dist.tsv"),
        overlap = os.path.join(STATS_DIR, "{sample}_overlap_stats.tsv")
    params:
        output_prefix = os.path.join(STATS_DIR, "{sample}")
    shell:
        """
        python scripts/gtf_stats.py \
            --reference {input.ref_gtf} \
            --query {input.query_gtf} \
            --output-prefix {params.output_prefix}
        """

