# File: scripts/calculate_pi.py

import argparse
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import glob
from collections import defaultdict

def calculate_pairwise_pi(seq1, seq2):
    """Calculate pairwise nucleotide diversity between two sequences"""
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

def calculate_between_groups_pi(sequences1, sequences2):
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

def read_genotype_matrix(matrix_file):
    """Read the genotype matrix to get allele counts for each ortholog"""
    df = pd.read_csv(matrix_file, sep="\t")
    
    # Define group1 (excluding RCC1749 and RCC3052)
    group1 = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC2482", 
              "RCC373", "RCC465", "RCC629", "RCC692", "RCC693", "RCC833"]
    
    # Count present (1) and absent (2) alleles for each ortholog in group1
    allele_counts = {}
    
    for ortholog_id in df['ortholog_id'].unique():
        ortholog_df = df[df['ortholog_id'] == ortholog_id]
        g1_rows = ortholog_df[ortholog_df['sample'].isin(group1)]
        
        # Count present alleles (presence = 1)
        present_count = sum(g1_rows['presence'] == 1)
        total_callable = sum(g1_rows['presence'] != 3)  # Exclude missing data
        
        if total_callable > 0:
            allele_freq = present_count / total_callable
            allele_counts[ortholog_id] = {
                'present_count': present_count,
                'total_callable': total_callable,
                'allele_freq': allele_freq
            }
    
    return allele_counts

def process_alignments(alignment_dir, genotype_matrix):
    """Process alignments to calculate Pi metrics for specific allele count groups"""
    # Read allele counts from genotype matrix
    allele_counts = read_genotype_matrix(genotype_matrix)
    
    # Group files by ortholog_id, side, and presence
    files = defaultdict(lambda: defaultdict(dict))
    
    # Find all separate alignment files for g1 group (Group1)
    for filepath in glob.glob(f"{alignment_dir}/separate/*.g1.*.mafft.fa"):
        filename = os.path.basename(filepath)
        parts = filename.split('.')
        ortholog_id = parts[0]
        side = parts[2].split('_')[0]  # left or right
        presence = parts[3]  # present or absent
        
        files[ortholog_id][side][presence] = filepath
    
    # Find all combined alignment files
    for filepath in glob.glob(f"{alignment_dir}/combined/*.combined.mafft.fa"):
        filename = os.path.basename(filepath)
        parts = filename.split('.')
        ortholog_id = parts[0]
        side = parts[1].split('_')[0]  # left or right
        
        files[ortholog_id][side]['combined'] = filepath
    
    # Results by ortholog and by count group
    all_results = []
    count_group_results = defaultdict(list)
    
    # Process each ortholog_id
    for ortholog_id in files:
        # Skip if no allele count information
        if ortholog_id not in allele_counts:
            continue
            
        present_count = allele_counts[ortholog_id]['present_count']
        total_callable = allele_counts[ortholog_id]['total_callable']
        
        # Skip singletons and fixed alleles
        if present_count == 1 or present_count == total_callable:
            continue
            
        # Keep only counts 2-10
        if present_count < 2 or present_count > 10:
            continue
            
        ortholog_results = {
            'ortholog_id': ortholog_id,
            'present_count': present_count,
            'total_callable': total_callable,
            'allele_freq': allele_counts[ortholog_id]['allele_freq']
        }
        
        # Calculate Pi for each side and take the average
        side_results = []
        
        for side in ['left', 'right']:
            if side not in files[ortholog_id]:
                continue
                
            side_data = files[ortholog_id][side]
            
            # Get sequences for present and absent groups
            present_seqs = []
            absent_seqs = []
            
            if 'present' in side_data:
                present_seqs = [str(record.seq) for record in 
                              SeqIO.parse(side_data['present'], 'fasta')]
            
            if 'absent' in side_data:
                absent_seqs = [str(record.seq) for record in 
                             SeqIO.parse(side_data['absent'], 'fasta')]
            
            # Calculate within-group Pi
            pi_present = calculate_group_pi(present_seqs)
            pi_absent = calculate_group_pi(absent_seqs)
            
            # Calculate between-group Pi from combined alignment
            pi_between = 0
            if 'combined' in side_data and present_seqs and absent_seqs:
                # Get record IDs for present sequences
                present_ids = set(r.id for r in SeqIO.parse(side_data['present'], 'fasta'))
                
                # Load combined alignment
                combined_records = list(SeqIO.parse(side_data['combined'], 'fasta'))
                
                # Separate into present and absent groups
                present_combined = [str(r.seq) for r in combined_records if r.id in present_ids]
                absent_combined = [str(r.seq) for r in combined_records if r.id not in present_ids]
                
                # Calculate between-group Pi
                pi_between = calculate_between_groups_pi(present_combined, absent_combined)
            
            side_results.append({
                'pi_present': pi_present,
                'pi_absent': pi_absent,
                'pi_between': pi_between
            })
        
        # Average results from both sides
        if side_results:
            for metric in ['pi_present', 'pi_absent', 'pi_between']:
                values = [result[metric] for result in side_results]
                ortholog_results[metric] = np.mean(values)
            
            # Add to results
            all_results.append(ortholog_results)
            
            # Add to count group results
            count_group_results[present_count].append(ortholog_results)
    
    # Compute summary statistics for each count group
    count_group_summary = []
    for count, orthologs in sorted(count_group_results.items()):
        if orthologs:
            summary = {
                'present_count': count,
                'loci_count': len(orthologs),
                'mean_pi_present': np.mean([o['pi_present'] for o in orthologs]),
                'mean_pi_absent': np.mean([o['pi_absent'] for o in orthologs]),
                'mean_pi_between': np.mean([o['pi_between'] for o in orthologs]),
                'median_pi_present': np.median([o['pi_present'] for o in orthologs]),
                'median_pi_absent': np.median([o['pi_absent'] for o in orthologs]),
                'median_pi_between': np.median([o['pi_between'] for o in orthologs]),
                'std_pi_present': np.std([o['pi_present'] for o in orthologs]),
                'std_pi_absent': np.std([o['pi_absent'] for o in orthologs]),
                'std_pi_between': np.std([o['pi_between'] for o in orthologs])
            }
            count_group_summary.append(summary)
    
    return pd.DataFrame(all_results), pd.DataFrame(count_group_summary)

def main():
    parser = argparse.ArgumentParser(description='Calculate Pi metrics for specific allele count groups (2-10)')
    parser.add_argument('--alignment_dir', required=True, help='Directory containing alignment files')
    parser.add_argument('--genotype_matrix', required=True, help='Genotype matrix file')
    parser.add_argument('--output_file', required=True, help='Output TSV file for Pi results')
    parser.add_argument('--summary_file', required=True, help='Output TSV file for count group summary')
    
    args = parser.parse_args()
    
    # Process alignments
    results_df, summary_df = process_alignments(args.alignment_dir, args.genotype_matrix)
    
    # Save results
    results_df.to_csv(args.output_file, sep='\t', index=False)
    print(f"Pi results saved to {args.output_file}")
    
    # Save count group summary
    summary_df.to_csv(args.summary_file, sep='\t', index=False)
    print(f"Count group summary saved to {args.summary_file}")

if __name__ == '__main__':
    main()