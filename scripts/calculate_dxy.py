# File: scripts/calculate_dxy.py

from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import numpy as np
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
    
    return presence_dict

def main():
    alignment_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    # Build presence dictionary
    presence_dict = build_presence_dict(alignment_dir)
    
    # Calculate dXY for each ortholog group
    results = []
    
    for ortholog_id in presence_dict:
        # Skip if missing data for either group
        if 'g1' not in presence_dict[ortholog_id] or 'g2' not in presence_dict[ortholog_id]:
            continue
            
        presence_g1 = presence_dict[ortholog_id]['g1']
        presence_g2 = presence_dict[ortholog_id]['g2']
        
        # Track metrics for each side (left/right)
        side_results = {}
        
        for side in ['left', 'right']:
            g1_file = f"{alignment_dir}/separate/{ortholog_id}.g1.{side}_flank.{presence_g1}.mafft.fa"
            g2_file = f"{alignment_dir}/separate/{ortholog_id}.g2.{side}_flank.{presence_g2}.mafft.fa"
            combined_file = f"{alignment_dir}/combined/{ortholog_id}.{side}_flank.combined.mafft.fa"
            
            # Skip if files don't exist
            if not os.path.exists(g1_file) or not os.path.exists(g2_file) or not os.path.exists(combined_file):
                continue
                
            # Load sequences
            g1_seqs = [str(record.seq) for record in SeqIO.parse(g1_file, 'fasta')]
            g2_seqs = [str(record.seq) for record in SeqIO.parse(g2_file, 'fasta')]
            
            # For between-group calculations, we need to identify sequences from combined alignment
            g1_ids = set(r.id for r in SeqIO.parse(g1_file, 'fasta'))
            combined_records = list(SeqIO.parse(combined_file, 'fasta'))
            
            g1_seqs_combined = [str(r.seq) for r in combined_records if r.id in g1_ids]
            g2_seqs_combined = [str(r.seq) for r in combined_records if r.id not in g1_ids]
            
            # Calculate nucleotide diversity metrics
            pi_g1 = calculate_group_pi(g1_seqs)
            pi_g2 = calculate_group_pi(g2_seqs)
            dxy_g1_vs_g2 = calculate_between_groups_dxy(g1_seqs_combined, g2_seqs_combined)
            
            side_results[side] = {
                'g1': pi_g1,
                'g2': pi_g2,
                'g1_vs_g2': dxy_g1_vs_g2
            }
        
        # Average results from both sides
        if side_results:
            result = {
                'ortholog_id': ortholog_id,
                'presence_g1': presence_g1,
                'presence_g2': presence_g2
            }
            
            for metric in ['g1', 'g2', 'g1_vs_g2']:
                values = [side_results.get(side, {}).get(metric, 0) for side in ['left', 'right'] 
                         if side in side_results]
                result[metric] = np.mean(values) if values else 0
                
            results.append(result)
    
    # Convert to DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(output_file, sep='\t', index=False)
    print(f"dXY results saved to {output_file}")

if __name__ == "__main__":
    main()