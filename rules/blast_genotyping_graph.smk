import os 

GENOME_DIR = "/scratch1/chris/introner-genotyping-pipeline/graph-assemblies/"
SEQUENCES = "/scratch1/chris/introner-genotyping-pipeline/introner_seq_truth_set/introner_truth_set.fa"
OUTPUT_DIR =  "/scratch1/chris/introner-genotyping-pipeline/graph_blast_results"
FINAL_SAMPLES = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
GTF_DIR = "/scratch1/chris/introner-genotyping-pipeline/gene_annotations/"
REF_DIR = "/scratch1/chris/introner-genotyping-pipeline/introner_files_for_Github/"

rule all:
    input:
        # expand("graph_blast_results/{genome_basename}_filtered_blast_results.unique.bed", genome_basename=SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.similarity_checked.fa"), genome_basename=FINAL_SAMPLES),
        expand(os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.similarity_checked.bed"), genome_basename=FINAL_SAMPLES)

# Rule to create a BLAST database for each genome
rule makeblastdb:
    input:
        genome= GENOME_DIR + "{genome_basename}.vg_paths.fa"
    output:
        db= GENOME_DIR + "{genome_basename}.nhr"  # A sample output from makeblastdb, there will be multiple files
    params:
        path= GENOME_DIR
    shell:
        """
        makeblastdb -in {input.genome} -dbtype nucl -out {params.path}{wildcards.genome_basename}
        """

# Rule to run BLASTn for each genome against the sequences
rule blast:
    input:
        db= GENOME_DIR + "{genome_basename}.nhr",  # Input BLAST database
        sequences=SEQUENCES  # Input sequences to search for
    output:
        os.path.join(OUTPUT_DIR, "{genome_basename}_blast_results.txt")
    params:
        db_prefix= GENOME_DIR + "{genome_basename}"  # Prefix of the BLAST database
    threads: 5
    shell:
        """
        blastn -num_threads {threads} -query {input.sequences} -db {params.db_prefix} -out {output} -outfmt "6 qseqid sseqid pident qlen slen length mismatch gapopen qstart qend sstart send evalue"
        """

# New rule to filter BLAST results for high confidence and length criteria
rule filter_blast_results:
    input:
        blast_results = os.path.join(OUTPUT_DIR, "{genome_basename}_blast_results.txt")
    output:
        high_conf_results = os.path.join(OUTPUT_DIR, "{genome_basename}_filtered_blast_results.txt")
    shell:
        """
        awk '($3 >= 80) && ($6 >= 0.8 * $4) && ($6 <= 1.2 * $4)' {input.blast_results} > {output.high_conf_results}
        """

rule convert_blast_to_bed:
    input:
        blast_results = os.path.join(OUTPUT_DIR, "{genome_basename}_filtered_blast_results.txt")
    output:
        bed = os.path.join(OUTPUT_DIR, "{genome_basename}_filtered_blast_results.bed")
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{start = ($11 < $12) ? $11 : $12; end = ($11 > $12) ? $11 : $12; print $2, start-1, end, $1, $3, $4, $5, $6, $7, $8, $9, $10, $13}}' {input.blast_results} > {output.bed}
        """

rule filter_bed_file:
    input:
        pre_bed = os.path.join(OUTPUT_DIR, "{genome_basename}_filtered_blast_results.bed")
    output:
        post_bed = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci.bed")
    shell:
        """
        python scripts/assess_blast_results_faster.py {input.pre_bed} {output.post_bed}
        """ 

rule make_flank_bed:
    input:
        bed = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci.bed")
    output:
        bed = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.bed")
    shell:
        """
        awk '{{print $1, $2 - 100, $3 + 100, $4}}' {input.bed} > {output.bed}
        """


rule bed_to_fa:
    input:
        bed = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.bed"),
        fa = os.path.join(GENOME_DIR, "{genome_basename}.vg_paths.fa")
    output:
        fa = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.fa"),
    shell:
        """
        awk 'OFS="\\t" {{print $1, $2, $3, $1":"$2"-"$3"|"$4}}' {input.bed} \
        | bedtools getfasta -fi {input.fa} -bed - -name \
        | sed '/^>/ s/::.*$//' > {output.fa}
        """

rule filter_repetitive_sequences:
    input:
        fa = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.fa")
    output:
        fa = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.fa")
    log: 
        os.path.join(OUTPUT_DIR, "logs", "{genome_basename}.removed_candidate_loci.log")
    shell:
        """
        python scripts/remove_low_complexity_sequences.py --input {input.fa} --output {output.fa} --log {log} --kmers 25
        """

rule check_sequence_similarity_to_reference:
    input:
        fa = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.fa"),
        gtf = os.path.join(GTF_DIR, "{genome_basename}.vg_paths.gtf"),
        ref_dir = REF_DIR,

    output:
        fa = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.similarity_checked.fa"),
        log1 = os.path.join(OUTPUT_DIR, "logs", "{genome_basename}.introner_similarity_check.log"),
        log2 = os.path.join(OUTPUT_DIR, "logs", "{genome_basename}.introner_similarity_check_summary.log")
    shell:
        """
        python scripts/validate_introners.py \
            --input_fasta {input.fa} \
            --gtf {input.gtf} \
            --ref_dir {input.ref_dir} \
            --output_fasta {output.fa} \
            --output_log {output.log1} \
            --similarity_cutoff 0.7
        """


rule make_filtered_bed:
    input:
        fa = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.similarity_checked.fa")
    output:
        bed = os.path.join(OUTPUT_DIR, "{genome_basename}.candidate_loci_plus_flanks.filtered.similarity_checked.bed")
    run:
        with open(input.fa, 'r') as f, open(output.bed, 'w') as o:
            for line in f:
                if line.startswith(">"):
                    try:
                        line = line[1:].strip()
                        region, metadata = line.split(" ")[0], line.split(" ")[1:]
                        region, identifier = region.split("|")
                        chrom, coords = region.split(":")
                        start, end = coords.split("-")
                        metadata = "\t".join(metadata)
                        o.write(f"{chrom}\t{start}\t{end}\t{identifier}\t{metadata}\n")
                    except ValueError as e:
                        raise ValueError(f"Error processing line: {line}. Metadata: {metadata}. Original error: {e}")

