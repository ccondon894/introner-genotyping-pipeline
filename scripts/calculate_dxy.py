from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import sys
import os
import glob

def calculate_pairwise_pi(seq1, seq2):
    """Calculate pairwise nucleotide diversity between two sequences, treating N's and gaps as missing data"""
    def is_valid_base(base):
        return base not in ['-', 'N', 'n']
    
    differences = sum(1 for a, b in zip(seq1, seq2) 
                     if a != b and is_valid_base(a) and is_valid_base(b))
    sites = sum(1 for a, b in zip(seq1, seq2) 
                if is_valid_base(a) and is_valid_base(b))
    return differences / sites if sites > 0 else 0

def calculate_group_pi(sequences):
    """Calculate average pi within a group of sequences"""
    if len(sequences) < 2:
        return 0
    
    n = len(sequences)
    total_pi = 0
    comparisons = 0
    
    for i in range(n):
        for j in range(i + 1, n):
            pi = calculate_pairwise_pi(sequences[i], sequences[j])
            total_pi += pi
            comparisons += 1
    
    return total_pi / comparisons if comparisons > 0 else 0

def calculate_between_groups_dxy(sequences1, sequences2):
    """Calculate average pi between two groups of sequences"""
    if not sequences1 or not sequences2:
        return 0
    
    total_pi = 0
    comparisons = 0
    
    for seq1 in sequences1:
        for seq2 in sequences2:
            pi = calculate_pairwise_pi(seq1, seq2)
            total_pi += pi
            comparisons += 1
    
    return total_pi / comparisons if comparisons > 0 else 0


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

mafft_directory = sys.argv[1]
output_file = sys.argv[2]

presence_dict = build_presence_dict(mafft_directory)

divergence_dict = {ortholog_id:  {'left': {'g1': 0, 'g2': 0, 'g1_vs_g2': 0},
                                'right': {'g1': 0, 'g2': 0, 'g1_vs_g2': 0}} for ortholog_id in presence_dict}

for ortholog_id in presence_dict:
    for side in ['left', 'right']:

        # make file paths
        g1_file = f"{mafft_directory}/separate/{ortholog_id}.g1.{side}_flank.{presence_dict[ortholog_id]['g1']}.mafft.fa"
        g2_file = f"{mafft_directory}/separate/{ortholog_id}.g2.{side}_flank.{presence_dict[ortholog_id]['g2']}.mafft.fa"
        g1_vs_g2_file = f"{mafft_directory}/combined/{ortholog_id}.{side}_flank.combined.mafft.fa"

        # load in aligned sequences
        g1_seqs = [str(record.seq) for record in SeqIO.parse(g1_file, 'fasta')]
        g2_seqs = [str(record.seq) for record in SeqIO.parse(g2_file, 'fasta')]
        g1_ids = set(r.id for r in SeqIO.parse(g1_file, 'fasta'))

        g1_seqs_combined = [str(r.seq) for r in SeqIO.parse(g1_vs_g2_file, 'fasta') if r.id in g1_ids]
        g2_seqs_combined = [str(r.seq) for r in SeqIO.parse(g1_vs_g2_file, 'fasta') if r.id not in g1_ids]

        # calculate pi and dxy
        pi_g1 = calculate_group_pi(g1_seqs)
        pi_g2 = calculate_group_pi(g2_seqs)
        dxy_g1_vs_g2 = calculate_between_groups_dxy(g1_seqs_combined, g2_seqs_combined)

        # record results
        divergence_dict[ortholog_id][side]['g1'] = pi_g1
        divergence_dict[ortholog_id][side]['g2'] = pi_g2
        divergence_dict[ortholog_id][side]['g1_vs_g2'] = dxy_g1_vs_g2

# average results from both sides
avg_divergence_dict = {ortholog_id : {'g1': 0,
                                      'g2': 0,
                                      'g1_vs_g2': 0}}
for ortholog_id in divergence_dict:
    avg_g1 = (divergence_dict[ortholog_id]['left']['g1'] + divergence_dict[ortholog_id]['right']['g1']) / 2
    avg_g2 = (divergence_dict[ortholog_id]['left']['g2'] + divergence_dict[ortholog_id]['right']['g2']) / 2
    avg_dxy = (divergence_dict[ortholog_id]['left']['g1_vs_g2'] + divergence_dict[ortholog_id]['right']['g1_vs_g2']) / 2
    presence_g1 = presence_dict[ortholog_id]['g1']
    presence_g2 = presence_dict[ortholog_id]['g2']
    avg_divergence_dict[ortholog_id] = {'g1': avg_g1, 'g2': avg_g2, 'g1_vs_g2': avg_dxy, 'presence_g1': presence_g1, 'presence_g2': presence_g2}


df = pd.DataFrame.from_dict(avg_divergence_dict, orient='index')
df['ortholog_id'] = df.index
df.reset_index(drop=True, inplace=True)
df = df[['ortholog_id', 'g1', 'g2', 'g1_vs_g2', 'presence_g1', 'presence_g2']]

df.to_csv(output_file, sep='\t', index=False)