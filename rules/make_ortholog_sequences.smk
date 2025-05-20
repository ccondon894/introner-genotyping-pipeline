OUTPUT_DIR = "/scratch1/chris/introner-genotyping-pipeline/flanking_pi_output"
ALIGNMENT_DIR = "/scratch1/chris/introner-genotyping-pipeline/mafft_flank_alignments"
GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"
ASSEMBLY_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies"
BAM_DIR = "/storage1/chris/introner-genotyping-bams"
DIVERSITY_DIR = "/scratch1/chris/introner-genotyping-pipeline/diversity_analysis"

import os
import glob
import json
from collections import defaultdict

all_samples = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
samples = [ "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
ref = ["CCMP1545"]
sides = ["left", "right"]

# Define the groups
group1 = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482", "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
group2 = ["RCC1749", "RCC3052"]

wildcard_constraints:
    asample = "|".join(all_samples),
    sample = "|".join(samples),
    ref = "|".join(ref),
    group = "g1|g2|g1_present|g1_absent|g1_combined"

# Function to get ortholog IDs after checkpoint execution
def get_ortholog_ids(wildcards=None):
    """Get list of ortholog IDs from processed files after checkpoint execution"""
    checkpoint_output = checkpoints.make_ortholog_fastas.get().output[0]
    
    # Look for actual files that were generated
    ortholog_ids = set()
    for fa_file in glob.glob(os.path.join(checkpoint_output, "*.fa")):
        # Extract ortholog ID from filename
        basename = os.path.basename(fa_file)
        ortholog_id = basename.split('.')[0]
        ortholog_ids.add(ortholog_id)
    
    return list(ortholog_ids)

# Function to get alignment files based on ortholog classification
def get_alignment_files(wildcards):
    """Get all alignment files that need to be created"""
    checkpoint_output = checkpoints.make_ortholog_fastas.get().output[0]
    
    # Load ortholog classification
    classification_file = os.path.join(ALIGNMENT_DIR, "ortholog_classification.json")
    with open(classification_file, "r") as f:
        classification = json.load(f)
    
    alignment_files = []
    
    # Process each ortholog based on its classification
    for oid, info in classification.items():
        category = info["category"]
        
        for side in sides:
            if category == "monomorphic_both_groups":
                # Both groups monomorphic - calculate Pi within each group and dXY between
                alignment_files.append(os.path.join(ALIGNMENT_DIR, f"{oid}.g1.{side}_flank.mafft.fa"))
                alignment_files.append(os.path.join(ALIGNMENT_DIR, f"{oid}.g2.{side}_flank.mafft.fa"))
                alignment_files.append(os.path.join(ALIGNMENT_DIR, f"{oid}.g1_g2_combined.{side}_flank.mafft.fa"))
            
            elif category == "polymorphic_group1":
                # Group1 polymorphic - calculate Pi within present, within absent, and between
                if info["g1_present_count"] > 1:
                    alignment_files.append(os.path.join(ALIGNMENT_DIR, f"{oid}.g1_present.{side}_flank.mafft.fa"))
                if info["g1_absent_count"] > 1:
                    alignment_files.append(os.path.join(ALIGNMENT_DIR, f"{oid}.g1_absent.{side}_flank.mafft.fa"))
                if info["g1_present_count"] > 0 and info["g1_absent_count"] > 0:
                    alignment_files.append(os.path.join(ALIGNMENT_DIR, f"{oid}.g1_combined.{side}_flank.mafft.fa"))
    
    return alignment_files

# Function to get all required files for the pipeline
def get_all_required_files(wildcards):
    """Get all files that need to be created by the pipeline"""
    files = []
    
    # Add consensus sequences
    files.extend(expand(os.path.join(OUTPUT_DIR, "{asample}.loci.{side}_flank.bed"), asample=all_samples, side=sides))
    files.extend(expand(os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.mpileup.bcf"), sample=samples, side=sides))
    files.extend(expand(os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.filtered.vcf.gz"), sample=samples, side=sides))
    files.extend(expand(os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.consensus.fa"), sample=samples, side=sides))
    files.extend(expand(os.path.join(OUTPUT_DIR, "{ref}.loci.{side}_flank.consensus.fa"), ref=ref, side=sides))
    
    # Add alignment files
    files.extend(get_alignment_files(wildcards))
    
    # Add diversity analysis output
    files.extend([
        os.path.join(DIVERSITY_DIR, "diversity_metrics.tsv"),
        os.path.join(DIVERSITY_DIR, "diversity_visualization.png")
    ])
    
    return files

rule all:
    input:
        get_all_required_files


rule make_bed_files:
    input:
        matrix = os.path.join(GT_DIR, "genotype_matrix_updated.tsv")
    output:
        left_bed = os.path.join(OUTPUT_DIR, "{asample}.loci.left_flank.bed"),
        right_bed = os.path.join(OUTPUT_DIR, "{asample}.loci.right_flank.bed")
    shell:
        """
        python scripts/make_sample_loci_flank_beds.py {input.matrix} {wildcards.asample} {output.left_bed} {output.right_bed}
        """

rule mpileup:
    input:
        bam = os.path.join(BAM_DIR, "{sample}.sorted.bam"),
        bed = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.bed"),
        fa = os.path.join(ASSEMBLY_DIR, "{sample}.vg_paths.fa")
    output:
        mpileup = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.mpileup.bcf")
    params:
        # Create a unique temp directory for this rule using wildcards
        temp_dir = os.path.join("temp_{sample}_{side}_sort")
    shell:
        """
        # Create temp directory
        mkdir -p {params.temp_dir}
        
        # Run commands
        bcftools mpileup \
            -d 100 -q 30 -Q 20 -A \
            -f {input.fa} -R {input.bed} \
            -Ou {input.bam} | \
        bcftools sort -T {params.temp_dir} -Ob -o {output.mpileup}
        
        # Cleanup temp directory
        rm -rf {params.temp_dir}
        """

rule call_variants:
    input:
        mpileup = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.mpileup.bcf"),
        fa = os.path.join(ASSEMBLY_DIR, "{sample}.vg_paths.fa")
    output:
        vcf = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.vcf.gz"),
        filt_vcf = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.filtered.vcf.gz")
    shell:
        """
        # Call variants
        bcftools call -m --ploidy 1 -Ou {input.mpileup} | \
        bcftools norm -f {input.fa} -m +both -Oz -o {output.vcf}

        # Filter variants
        bcftools filter -i 'DP>=10' {output.vcf} -o {output.filt_vcf}
        tabix -p vcf {output.filt_vcf}
        """

rule build_sample_consensus:
    input:
        filt_vcf = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.filtered.vcf.gz"),
        bed = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.bed"),
        fa = os.path.join(ASSEMBLY_DIR, "{sample}.vg_paths.fa")
    output:
        cons = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.consensus.fa")
    params:
        cons_genome = os.path.join(OUTPUT_DIR, "{sample}.{side}.consensus.fa"),
        loci = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.txt")
    shell:
        """
        # Generate consensus genome
        bcftools consensus -a N -f {input.fa} -o {params.cons_genome} {input.filt_vcf}

        # Index the consensus genome
        samtools faidx {params.cons_genome}

        # Extract loci from BED and generate a loci list
        awk '{{print $1":"$2"-"$3}}' {input.bed} > {params.loci}

        # Extract the consensus sequences for the given loci
        samtools faidx {params.cons_genome} $(cat {params.loci}) > {output.cons}

        # Clean up temporary files
        rm -f {params.cons_genome} {params.cons_genome}.fai {params.loci}
        """

rule build_ref_consensus:
    input:
        bed = os.path.join(OUTPUT_DIR, "{ref}.loci.{side}_flank.bed"),
        fa = os.path.join(ASSEMBLY_DIR, "{ref}.vg_paths.fa")
    output:
        fa = os.path.join(OUTPUT_DIR, "{ref}.loci.{side}_flank.consensus.fa")
    shell:
        """
        bedtools getfasta -fi {input.fa} -bed {input.bed} -fo {output.fa} -name 
        sed -i 's/ortholog_id_[0-9]\{{4\}}:://g' {output.fa}
        """

checkpoint make_ortholog_fastas:
    input:
        fastas = expand(os.path.join(OUTPUT_DIR, "{asample}.loci.{side}_flank.consensus.fa"), asample=all_samples, side=sides),
        gt_matrix = os.path.join(GT_DIR, "genotype_matrix_updated.tsv"),
        fasta_dir = OUTPUT_DIR
    output:
        alignment_dir = directory(ALIGNMENT_DIR),
        classification = os.path.join(ALIGNMENT_DIR, "ortholog_classification.json")
    params:
        group1_str = ",".join(group1),
        group2_str = ",".join(group2)
    shell:
        """
        mkdir -p {output.alignment_dir}
        python scripts/classify_and_prepare_orthologs.py \
            {input.gt_matrix} {input.fasta_dir} {output.alignment_dir} \
            --group1 {params.group1_str} --group2 {params.group2_str}
        """

rule align_monomorphic_group:
    input:
        fasta = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.{group}.{side}_flank.fa")
    output:
        aligned = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.{group}.{side}_flank.mafft.fa")
    wildcard_constraints:
        group = "g1|g2|g1_present|g1_absent"
    shell:
        """
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}
        """

rule align_combined_groups:
    input:
        g1 = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1.{side}_flank.fa"),
        g2 = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g2.{side}_flank.fa")
    output:
        combined = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1_g2_combined.{side}_flank.fa"),
        aligned = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1_g2_combined.{side}_flank.mafft.fa")
    shell:
        """
        # Combine the sequences
        cat {input.g1} {input.g2} > {output.combined}
        
        # Align the combined sequences
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {output.combined} > {output.aligned}
        """

rule align_combined_polymorphic:
    input:
        present = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1_present.{side}_flank.fa"),
        absent = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1_absent.{side}_flank.fa")
    output:
        combined = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1_combined.{side}_flank.fa"),
        aligned = os.path.join(ALIGNMENT_DIR, "{ortholog_id}.g1_combined.{side}_flank.mafft.fa")
    shell:
        """
        # Combine the sequences
        cat {input.present} {input.absent} > {output.combined}
        
        # Align the combined sequences
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {output.combined} > {output.aligned}
        """

rule calculate_diversity_metrics:
    input:
        alignment_files = get_alignment_files
    output:
        metrics = os.path.join(DIVERSITY_DIR, "diversity_metrics.tsv")
    params:
        alignment_dir = ALIGNMENT_DIR,
        classification = os.path.join(ALIGNMENT_DIR, "ortholog_classification.json")
    shell:
        """
        mkdir -p {DIVERSITY_DIR}
        python scripts/calculate_diversity_metrics.py \
            --alignment_dir {params.alignment_dir} \
            --classification {params.classification} \
            --output {output.metrics}
        """

rule visualize_diversity:
    input:
        metrics = os.path.join(DIVERSITY_DIR, "diversity_metrics.tsv")
    output:
        plot = os.path.join(DIVERSITY_DIR, "diversity_visualization.png")
    shell:
        """
        python scripts/visualize_diversity.py \
            --input {input.metrics} \
            --output {output.plot}
        """



