"""
Snakemake workflow for genome comparison and variant calling.
This workflow:
1. Aligns a query genome to a reference genome using Minimap2
2. Calls variants using bcftools
3. Counts SNPs and indels in the resulting VCF file
"""

# Configuration
configfile: "config.yaml"

# Default values if config is not provided
REF_GENOME = config.get("reference_genome", "reference.fa")
QUERY_GENOME = config.get("query_genome", "query.fa")
DIVERGENCE = config.get("expected_divergence", "5")  # Default to asm5 preset
OUTPUT_DIR = config.get("output_dir", "results")

# Final output target
rule all:
    input:
        f"{OUTPUT_DIR}/variant_counts.txt"

# Align query genome to reference using Minimap2
rule align_genomes:
    input:
        ref = REF_GENOME,
        query = QUERY_GENOME
    output:
        sam = temp(f"{OUTPUT_DIR}/alignment.sam")
    params:
        preset = f"asm{DIVERGENCE}"
    threads: config.get("threads", 8)
    log:
        f"{OUTPUT_DIR}/logs/minimap2.log"
    shell:
        """
        minimap2 -ax {params.preset} -t {threads} {input.ref} {input.query} > {output.sam} 2> {log}
        """

# Convert SAM to sorted BAM
rule sam_to_bam:
    input:
        sam = f"{OUTPUT_DIR}/alignment.sam"
    output:
        bam = f"{OUTPUT_DIR}/alignment.sorted.bam",
        bai = f"{OUTPUT_DIR}/alignment.sorted.bam.bai"
    threads: config.get("threads", 8)
    log:
        f"{OUTPUT_DIR}/logs/samtools.log"
    shell:
        """
        samtools view -bS {input.sam} | samtools sort -@ {threads} -o {output.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """

# Call variants using bcftools
rule call_variants:
    input:
        ref = REF_GENOME,
        bam = f"{OUTPUT_DIR}/alignment.sorted.bam"
    output:
        vcf = f"{OUTPUT_DIR}/variants.vcf"
    log:
        f"{OUTPUT_DIR}/logs/bcftools.log"
    shell:
        """
        bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv -Ov -o {output.vcf} 2> {log}
        """

# Count SNPs and indels in the VCF file
rule count_variants:
    input:
        vcf = f"{OUTPUT_DIR}/variants.vcf"
    output:
        counts = f"{OUTPUT_DIR}/variant_counts.txt"
    log:
        f"{OUTPUT_DIR}/logs/count_variants.log"
    shell:
        """
        echo "Genome Comparison Results" > {output.counts}
        echo "------------------------" >> {output.counts}
        echo "Reference genome: {REF_GENOME}" >> {output.counts}
        echo "Query genome: {QUERY_GENOME}" >> {output.counts}
        echo "------------------------" >> {output.counts}
        echo -n "SNPs: " >> {output.counts}
        grep -v "^#" {input.vcf} | awk '$5 != "." && length($4) == 1 && length($5) == 1' | wc -l >> {output.counts}
        
        echo -n "Insertions: " >> {output.counts}
        grep -v "^#" {input.vcf} | awk '$5 != "." && length($4) == 1 && length($5) > 1' | wc -l >> {output.counts}
        
        echo -n "Deletions: " >> {output.counts}
        grep -v "^#" {input.vcf} | awk '$5 != "." && length($4) > 1 && length($5) == 1' | wc -l >> {output.counts}
        
        # Calculate total variants
        echo -n "Total variants: " >> {output.counts}
        grep -v "^#" {input.vcf} | wc -l >> {output.counts}
        
        # Optional: Add alignment stats
        echo "------------------------" >> {output.counts}
        echo "Alignment statistics:" >> {output.counts}
        samtools flagstat {OUTPUT_DIR}/alignment.sorted.bam >> {output.counts} 2> {log}
        """