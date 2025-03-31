import glob
import os

# Config
FLANK_SIDES = ['left', 'right']
PRESENCE_STATES = ['present', 'absent']
SEQ_DIR = "/scratch1/chris/introner-genotyping-pipeline/mafft_flank_alignments"
WORKING_DIR = "scratch1/chris/introner-genotyping-pipeline"

# Get all ortholog IDs from the input files
def get_ortholog_ids():
    files = glob.glob(os.path.join(SEQ_DIR, "*.fa"))
    return list(set([f.split("/")[-1].split('.')[0] for f in files]))

ORTHOLOG_IDS = get_ortholog_ids()

rule all:
    input:
        os.path.join(SEQ_DIR, "diversity_metrics.tsv")

rule combine_sequences:
    input:
        present = os.path.join(SEQ_DIR, "{ortholog_id}.{side}_flank.present.fa"),
        absent = os.path.join(SEQ_DIR, "{ortholog_id}.{side}_flank.absent.fa")
    output:
        combined = os.path.join(SEQ_DIR, "combined/{ortholog_id}.{side}_flank.combined.fa")
    shell:
        "cat {input.present} {input.absent} > {output.combined}"

rule align_separate:
    input:
        fasta = os.path.join(SEQ_DIR, "{ortholog_id}.{side}_flank.{presence}.fa")
    output:
        aligned = os.path.join(SEQ_DIR, "separate/{ortholog_id}.{side}_flank.{presence}.mafft.fa")
    shell:
        "mafft {input.fasta} > {output.aligned}"

rule align_combined:
    input:
        fasta = os.path.join(SEQ_DIR, "combined/{ortholog_id}.{side}_flank.combined.fa")
    output:
        aligned = os.path.join(SEQ_DIR, "combined/{ortholog_id}.{side}_flank.combined.mafft.fa")
    shell:
        "mafft --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}"

rule calculate_pi:
    input:
        separate_alignments = expand(
            os.path.join(SEQ_DIR, "separate/{ortholog_id}.{side}_flank.{presence}.mafft.fa"),
            ortholog_id=ORTHOLOG_IDS,
            side=FLANK_SIDES,
            presence=PRESENCE_STATES
        ),
        combined_alignments = expand(
            os.path.join(SEQ_DIR, "combined/{ortholog_id}.{side}_flank.combined.mafft.fa"),
            ortholog_id=ORTHOLOG_IDS,
            side=FLANK_SIDES
        )
    output:
        metrics = os.path.join(SEQ_DIR, "diversity_metrics.tsv")
    shell:
        """
        echo '{input.separate_alignments} {input.combined_alignments}' | tr ' ' '\n' > input_files.txt && \
        python scripts/calculate_pi.py --input_list input_files.txt --output_file {output.metrics}
        """