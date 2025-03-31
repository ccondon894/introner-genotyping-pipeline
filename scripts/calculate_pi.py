from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from collections import defaultdict
import pandas as pd
from pathlib import Path
import argparse

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

def process_alignments(input_files):
    results = []
    
    # Group files by ortholog_id and side
    alignments = defaultdict(lambda: defaultdict(dict))
    for filepath in input_files:
        filename = Path(filepath).name
        parts = filename.split('.')
        ortholog_id = parts[0]
        side = parts[1].split('_')[0]
        if 'combined' in filename:
            alignments[ortholog_id][side]['combined'] = filepath
        else:
            presence = parts[2]
            alignments[ortholog_id][side][presence] = filepath

    # Process each ortholog
    for ortholog_id in alignments:
        ortholog_results = {'ortholog_id': ortholog_id}
        
        # Calculate pi for each side and take the average
        side_results = []
        for side in ['left', 'right']:
            if side not in alignments[ortholog_id]:
                continue
                
            # Calculate within-group pi from separate alignments
            present_seqs = []
            absent_seqs = []
            
            if 'present' in alignments[ortholog_id][side]:
                present_seqs = [str(record.seq) for record in 
                              SeqIO.parse(alignments[ortholog_id][side]['present'], 'fasta')]
            
            if 'absent' in alignments[ortholog_id][side]:
                absent_seqs = [str(record.seq) for record in 
                             SeqIO.parse(alignments[ortholog_id][side]['absent'], 'fasta')]
            
            pi_present = calculate_group_pi(present_seqs)
            pi_absent = calculate_group_pi(absent_seqs)
            
            # Calculate between-group pi from combined alignment
            if 'combined' in alignments[ortholog_id][side]:
                records = list(SeqIO.parse(alignments[ortholog_id][side]['combined'], 'fasta'))
                # Identify present and absent sequences in combined alignment
                # Assuming sequence headers are consistent between files
                present_ids = set(r.id for r in SeqIO.parse(alignments[ortholog_id][side]['present'], 'fasta'))
                present_seqs_combined = [str(r.seq) for r in records if r.id in present_ids]
                absent_seqs_combined = [str(r.seq) for r in records if r.id not in present_ids]
                pi_between = calculate_between_groups_pi(present_seqs_combined, absent_seqs_combined)
            else:
                pi_between = 0
            
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
        
        results.append(ortholog_results)
    
    return pd.DataFrame(results)

def main():
    parser = argparse.ArgumentParser(description='Calculate pi metrics for sequence alignments')
    parser.add_argument('--input_list', required=True,
                      help='File containing list of aligned fasta files')
    parser.add_argument('--output_file', required=True,
                      help='Output TSV file')
    
    args = parser.parse_args()
    
    # Read input files from list
    with open(args.input_list) as f:
        input_files = [line.strip() for line in f]
    
    # Process alignments
    df = process_alignments(input_files)
    
    # Ensure columns are in the correct order
    df = df[['ortholog_id', 'pi_present', 'pi_absent', 'pi_between']]
    
    # Save results
    df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()