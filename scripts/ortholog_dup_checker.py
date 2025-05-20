#!/usr/bin/env python3

import pandas as pd
import argparse
from collections import defaultdict

def check_unique_sequences(input_file):
    """
    Check if each sequence in the genotype matrix appears only once.
    Sequences are identified by their sample, start, end, and contig values.
    Rows with presence = 3 (missing data) are excluded.
    
    Args:
        input_file: Path to the genotype matrix file
    
    Returns:
        A tuple of (is_unique, duplicates) where duplicates is a dict of duplicate sequences
    """
    # Read the genotype matrix
    df = pd.read_csv(input_file, sep='\t')
    
    # Filter out rows with presence = 3
    df = df[df['presence'] != 3]
    
    # Create a dictionary to track occurrences
    sequence_counts = defaultdict(list)
    
    # Count occurrences of each sequence
    for idx, row in df.iterrows():
        # Create a tuple that uniquely identifies the sequence
        seq_key = (row['sample'], row['start'], row['end'], row['contig'])
        sequence_counts[seq_key].append(row['ortholog_id'])
    
    # Find duplicates (sequences with more than one occurrence)
    duplicates = {seq: orthologs for seq, orthologs in sequence_counts.items() if len(orthologs) > 1}
    
    # Return results
    return len(duplicates) == 0, duplicates

def main():
    parser = argparse.ArgumentParser(description='Check for duplicate sequences in genotype matrix')
    parser.add_argument('input_file', help='Path to input genotype matrix TSV file')
    parser.add_argument('--output', help='Optional output file for detailed results')
    
    args = parser.parse_args()
    
    # Check for duplicates
    is_unique, duplicates = check_unique_sequences(args.input_file)
    
    # Print results to console
    if is_unique:
        print("SUCCESS: All sequences in the genotype matrix are unique!")
    else:
        print(f"WARNING: Found {len(duplicates)} sequences that appear multiple times.")
        
    # Write detailed results to output file if specified
    if args.output and duplicates:
        with open(args.output, 'w') as f:
            f.write("sample\tstart\tend\tcontig\tortholog_ids\n")
            for (sample, start, end, contig), orthologs in duplicates.items():
                f.write(f"{sample}\t{start}\t{end}\t{contig}\t{','.join(orthologs)}\n")
        print(f"Detailed results written to {args.output}")

if __name__ == "__main__":
    main()