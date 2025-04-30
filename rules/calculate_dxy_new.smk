import os
import glob
import json

SEQ_DIR = "/scratch1/chris/introner-genotyping-pipeline/flanking_pi_output"
ALIGNMENT_DIR = "/scratch1/chris/introner-genotyping-pipeline/mafft_flank_alignments_dxy"
GT_DIR = "/scratch1/chris/introner-genotyping-pipeline/genotype_matrixes"

FLANK_SIDES = ['left', 'right']
PRESENCE_STATES = ['present', 'absent']
GROUPS = ['g1', 'g2']

# Function to build presence dict
def build_presence_dict(output_dir):
    """Build dictionary of presence/absence for g1/g2 pairs by examining files in directory"""
    presence_dict = {}
    
    # Get all relevant files
    files = glob.glob(f"{output_dir}/*.fa")
    
    for f in files:
        # Parse filename to get components
        basename = os.path.basename(f)
        parts = basename.split('.')
        if len(parts) >= 4:  # Make sure we have enough parts
            ortholog_id = parts[0]  # First element is ortholog_id
            group = parts[1]        # g1 or g2
            presence = parts[3]     # Get presence value
            
            # Initialize nested dict if ortholog_id not seen
            if ortholog_id not in presence_dict:
                presence_dict[ortholog_id] = {}
                
            # Add presence value for this group
            presence_dict[ortholog_id][group] = presence
    
    # Save the dict for future reference
    with open(os.path.join(output_dir, "presence_dict.json"), "w") as f:
        json.dump(presence_dict, f)
    
    return presence_dict

# Function to get ortholog IDs after checkpoint
def get_ortholog_ids(wildcards=None):
    checkpoint_output = checkpoints.make_dxy_fastas.get().output[0]
    
    # Load presence dict that was created by the checkpoint
    presence_dict_file = os.path.join(ALIGNMENT_DIR, "presence_dict.json")
    if os.path.exists(presence_dict_file):
        with open(presence_dict_file, "r") as f:
            presence_dict = json.load(f)
        return list(presence_dict.keys())
    else:
        # Fallback if file doesn't exist (shouldn't happen in normal flow)
        files = glob.glob(os.path.join(ALIGNMENT_DIR, "*.g1.*.fa"))
        return list(set([os.path.basename(f).split('.')[0] for f in files]))

# Function to get sequence inputs based on presence dict
def get_sequence_input(wildcards):
    # Load presence dict
    presence_dict_file = os.path.join(ALIGNMENT_DIR, "presence_dict.json")
    with open(presence_dict_file, "r") as f:
        presence_dict = json.load(f)
    
    return {
        "g1": f"{ALIGNMENT_DIR}/{wildcards.ortholog_id}.g1.{wildcards.side}_flank.{presence_dict[wildcards.ortholog_id]['g1']}.fa",
        "g2": f"{ALIGNMENT_DIR}/{wildcards.ortholog_id}.g2.{wildcards.side}_flank.{presence_dict[wildcards.ortholog_id]['g2']}.fa"
    }

# Function to get input for separate alignment
def get_separate_input(wildcards):
    presence_dict_file = os.path.join(ALIGNMENT_DIR, "presence_dict.json")
    with open(presence_dict_file, "r") as f:
        presence_dict = json.load(f)
    
    # The filename pattern should be adjusted based on your actual files
    return f"{ALIGNMENT_DIR}/{wildcards.ortholog_id}.{wildcards.group}.{wildcards.side}_flank.{presence_dict[wildcards.ortholog_id][wildcards.group]}.fa"

# Function to get all MAFFT files after they're created
def get_all_mafft_files(wildcards):
    checkpoint_output = checkpoints.make_dxy_fastas.get().output[0]
    
    # Get ortholog IDs
    ortholog_ids = get_ortholog_ids()
    
    # Generate a list of expected files
    mafft_files = []
    presence_dict_file = os.path.join(ALIGNMENT_DIR, "presence_dict.json")
    with open(presence_dict_file, "r") as f:
        presence_dict = json.load(f)
    
    # Add combined alignments
    for oid in ortholog_ids:
        for side in FLANK_SIDES:
            mafft_files.append(os.path.join(ALIGNMENT_DIR, f"combined/{oid}.{side}_flank.combined.mafft.fa"))
    
    # Add separate alignments
    for oid in ortholog_ids:
        for group in GROUPS:
            for side in FLANK_SIDES:
                presence = presence_dict[oid][group]
                mafft_files.append(os.path.join(ALIGNMENT_DIR, f"separate/{oid}.{group}.{side}_flank.{presence}.mafft.fa"))
    
    return mafft_files

# Main workflow
rule all:
    input:
        # Final output - changed to be outside the directory
        os.path.join(os.path.dirname(ALIGNMENT_DIR), "dxy_metrics.tsv")

# First rule to create FASTA files and directory structure
checkpoint make_dxy_fastas:
    input:
        metadata = os.path.join(GT_DIR, "genotype_matrix.metadata.4.updated.tsv"),
        fasta_dir = SEQ_DIR,
        calls = os.path.join(GT_DIR, "genotype_matrix.calls.4.updated.tsv")
    output:
        data_dir = directory(os.path.join(ALIGNMENT_DIR, "data")),  # Changed output to a subdirectory
        flag = os.path.join(ALIGNMENT_DIR, "setup_complete.flag")   # Added a flag file as output
    run:
        # Create main directory
        if not os.path.exists(ALIGNMENT_DIR):
            shell("mkdir -p {ALIGNMENT_DIR}")
        
        # Create data subdirectory
        if not os.path.exists(output.data_dir):
            shell("mkdir -p {output.data_dir}")
        
        # Create other subdirectories
        shell("mkdir -p {ALIGNMENT_DIR}/combined {ALIGNMENT_DIR}/separate")
        
        # Run script to create fastas
        shell("python scripts/make_ortholog_group_fastas_dxy.py {input.metadata} {input.calls} {input.fasta_dir} {ALIGNMENT_DIR}")
        
        # Build and save presence dict
        build_presence_dict(ALIGNMENT_DIR)
        
        # Create flag file
        shell("touch {output.flag}")

rule combine_sequences:
    input:
        flag = os.path.join(ALIGNMENT_DIR, "setup_complete.flag"),
        g1 = lambda wildcards: get_sequence_input(wildcards)["g1"],
        g2 = lambda wildcards: get_sequence_input(wildcards)["g2"]
    output:
        combined = os.path.join(ALIGNMENT_DIR, "combined/{ortholog_id}.{side}_flank.combined.fa")
    shell:
        """
        mkdir -p $(dirname {output.combined})
        cat {input.g1} {input.g2} > {output.combined}
        """

# Rule to align separate sequences
rule align_separate:
    input:
        flag = os.path.join(ALIGNMENT_DIR, "setup_complete.flag"),
        fasta = get_separate_input
    output:
        aligned = os.path.join(ALIGNMENT_DIR, "separate/{ortholog_id}.{group}.{side}_flank.{presence}.mafft.fa")
    shell:
        """
        mkdir -p $(dirname {output.aligned})
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}
        """

# Rule to align combined sequences
rule align_combined:
    input:
        flag = os.path.join(ALIGNMENT_DIR, "setup_complete.flag"),
        fasta = os.path.join(ALIGNMENT_DIR, "combined/{ortholog_id}.{side}_flank.combined.fa")
    output:
        aligned = os.path.join(ALIGNMENT_DIR, "combined/{ortholog_id}.{side}_flank.combined.mafft.fa")
    shell:
        """
        mkdir -p $(dirname {output.aligned})
        mafft --adjustdirection --maxiterate 1000 --globalpair --quiet {input.fasta} > {output.aligned}
        """

# Rule to calculate Dxy
rule calculate_dxy:
    input:
        flag = os.path.join(ALIGNMENT_DIR, "setup_complete.flag"),
        mafft_files = get_all_mafft_files
    params:
        alignment_dir = ALIGNMENT_DIR
    output:
        # Output file is now outside the ALIGNMENT_DIR
        metrics = os.path.join(os.path.dirname(ALIGNMENT_DIR), "dxy_metrics.tsv")
    shell:
        """
        python scripts/calculate_dxy.py {params.alignment_dir} {output.metrics}
        """