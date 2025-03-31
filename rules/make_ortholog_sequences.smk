
OUTPUT_DIR = "/scratch1/chris/introner-genotyping-pipeline/flanking_pi_output"
ALIGNMENT_DIR = "/scratch1/chris/introner-genotyping-pipeline/mafft_flank_alignments"
GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"
ASSEMBLY_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies"
BAM_DIR = "/storage1/chris/introner-genotyping"

all_samples = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]
samples = ["CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]
ref = ["CCMP1545"]
sides = ["left", "right"]

wildcard_constraints:
    asample = "|".join(all_samples),
    sample = "|".join(samples),
    ref = "|".join(ref)


def get_ortholog_files(wildcards):
    checkpoint_output = checkpoints.make_ortholog_fastas.get().output[0]
    
    # Instead of using glob_wildcards, let's directly get the files
    import glob
    import os
    
    # Get actual files that exist
    existing_files = glob.glob(os.path.join(checkpoint_output, "*.fa"))
    
    # Create the corresponding mafft output files
    mafft_files = [f.replace(".fa", ".mafft.fa") for f in existing_files]
    
    print(f"Found {len(existing_files)} input files")
    print(f"Requesting {len(mafft_files)} MAFFT alignments")
    
    return mafft_files


rule all:
    input:
        expand(os.path.join(OUTPUT_DIR, "{asample}.loci.{side}_flank.bed"), asample=all_samples, side=sides),
        expand(os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.mpileup.bcf"), sample=samples, side=sides),
        expand(os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.filtered.vcf.gz"),sample=samples, side=sides),
        expand(os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.consensus.fa"), sample=samples, side=sides),
        expand(os.path.join(OUTPUT_DIR, "{ref}.loci.{side}_flank.consensus.fa"), ref=ref, side=sides),
        ALIGNMENT_DIR


rule make_bed_files:
    input:
        matrix = os.path.join(GT_DIR, "genotype_matrix.metadata.3.updated.tsv")
    output:
        left_bed = os.path.join(OUTPUT_DIR, "{asample}.loci.left_flank.bed"),
        right_bed = os.path.join(OUTPUT_DIR, "{asample}.loci.right_flank.bed")
    shell:
        """
        python scripts/make_sample_loci_flank_beds.py {input.matrix} {wildcards.asample} {output.left_bed} {output.right_bed}
        """

rule mpileup:
    input:
        bam = os.path.join(BAM_DIR, "{sample}.surject_to.{sample}.paired.filtered.bam"),
        bed = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.bed"),
        fa = os.path.join(ASSEMBLY_DIR, "{sample}.vg_paths.fa")
    output:
        mpileup = os.path.join(OUTPUT_DIR, "{sample}.loci.{side}_flank.mpileup.bcf")
    shell:
        """
        bcftools mpileup \
            -d 100 -q 30 -Q 20 -A \
            -f {input.fa} -R {input.bed} \
            -Ou {input.bam} | \
        bcftools sort -T temp -Ob -o {output.mpileup}
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
       metadata = os.path.join(GT_DIR, "genotype_matrix.metadata.3.updated.tsv"),
       fasta_dir = OUTPUT_DIR,
       calls = os.path.join(GT_DIR, "genotype_matrix.calls.3.updated.tsv")
   params:
       fasta_dir = ALIGNMENT_DIR
   output:
       directory(ALIGNMENT_DIR)
   shell:
       """
       mkdir -p {output}
       python scripts/make_ortholog_group_fastas.py {input.metadata} {input.calls} {input.fasta_dir} {params.fasta_dir}
       """



