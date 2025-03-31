import os
import glob

SEQ_DIR = "/scratch1/chris/introner-genotyping-pipeline/flanking_pi_output"
ALIGNMENT_DIR = "/scratch1/chris/introner-genotyping-pipeline/mafft_flank_alignments_dxy"
GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"

FLANK_SIDES = ['left', 'right']
PRESENCE_STATES = ['present', 'absent']
GROUPS = ['g1', 'g2']
# Get all ortholog IDs from the input files
def get_ortholog_ids():
    checkpoint_output = checkpoints.build_presence_dict.get().output[0]
    files = glob.glob(os.path.join(ALIGNMENT_DIR, "*.fa"))
    return list(set([f.split("/")[-1].split('.')[0] for f in files]))

def build_presence_dict(output_dir):
    """Build dictionary of presence/absence for g1/g2 pairs by examining files in directory"""
    
    presence_dict = {}
    
    # Get all relevant files
    files = glob.glob(f"{output_dir}/*.fa")
    
    for f in files:
        # Parse filename to get components
        basename = os.path.basename(f)
        ortholog_id = basename.split('.')[0]  # Assuming first element is ortholog_id
        group = basename.split('.')[1]        # g1 or g2
        presence = basename.split('.')[3]  # Get presence value
        
        # Initialize nested dict if ortholog_id not seen
        if ortholog_id not in presence_dict:
            presence_dict[ortholog_id] = {}
            
        # Add presence value for this group
        presence_dict[ortholog_id][group] = presence

    return presence_dict

def get_mafft_files(output_dir):

    separate_files = glob.glob(f"{output_dir}/separate/*.mafft.fa")
    combined_files = glob.glob(f"{output_dir}/combined/*.mafft.fa")
    all_files = separate_files + combined_files
    return all_files


rule all:
    input:
        ALIGNMENT_DIR,
        expand(os.path.join(ALIGNMENT_DIR, "combined/{ortholog_id}.{side}_flank.combined.fa"), ortholog_id=ORTHOLOG_IDS, side=FLANK_SIDES),
        

rule make_dxy_fastas:
    input:
        metadata = os.path.join(GT_DIR, "genotype_matrix.metadata.3.updated.tsv"),
        fasta_dir = SEQ_DIR,
        calls = os.path.join(GT_DIR, "genotype_matrix.calls.3.updated.tsv")
    params:
        fasta_dir = ALIGNMENT_DIR
    output:
        directory(ALIGNMENT_DIR)
    shell:
        """
        mkdir -p {output}
        python scripts/make_ortholog_group_fastas_dxy.py {input.metadata} {input.calls} {input.fasta_dir} {params.fasta_dir}
        """


checkpoint build_presence_dict:
    input:
        alignment_dir = ALIGNMENT_DIR
    output:
        temp(os.path.join(ALIGNMENT_DIR, "presence_dict.done"))
    run:
        presence_dict = build_presence_dict(input.alignment_dir)
        with open(output[0], "w") as f:
            f.write("Done")

ORTHOLOG_IDS = get_ortholog_ids()
print(presence_dict)



rule combine_sequences:
    input:
        g1 = lambda wildcards: f"{ALIGNMENT_DIR}/{wildcards.ortholog_id}.g1.{wildcards.side}_flank.{presence_dict[wildcards.ortholog_id]['g1']}.fa",
        g2 = lambda wildcards: f"{ALIGNMENT_DIR}/{wildcards.ortholog_id}.g2.{wildcards.side}_flank.{presence_dict[wildcards.ortholog_id]['g2']}.fa"
    output:
        combined = os.path.join(ALIGNMENT_DIR, "combined/{ortholog_id}.{side}_flank.combined.fa")
    shell:
        """
        cat {input.g1} {input.g2} > {output.combined}
        """


rule align_separate:
    input:
        fasta = lambda wildcards: f"{ALIGNMENT_DIR}/separate/{wildcards.ortholog_id}.{wildcards.group}.{wildcards.side}_flank.{presence_dict[wildcards.ortholog_id][wildcards.group]}.fa",
    params:
        presence = lambda wildcards: presence_dict[wildcards.ortholog_id][wildcards.group]
    output:
        aligned = os.path.join(ALIGNMENT_DIR, "separate/{ortholog_id}.{group}.{side}_flank.{presence}.mafft.fa")
    shell:
        "mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}"

rule align_combined:
    input:
        fasta = lambda wildcards: f"{ALIGNMENT_DIR}/combined/{wildcards.ortholog_id}.{wildcards.side}_flank.combined.fa",
    output:
        aligned = os.path.join(ALIGNMENT_DIR, "combined/{ortholog_id}.{side}_flank.combined.mafft.fa")
    shell:
        "mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}"


all_files = get_mafft_files(ALIGNMENT_DIR)

rule calculate_dxy:
    input:
        mafft_files = all_files
    output:
        metrics = f"{ALIGNMENT_DIR}/dxy_metrics.tsv"
    shell:
        """
        python scripts/calculate_dxy.py
        """

