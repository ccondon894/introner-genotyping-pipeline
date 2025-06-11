#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd
import json
from Bio import SeqIO
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Classify and prepare ortholog groups for diversity analysis')
    parser.add_argument('genotype_matrix', help='Path to genotype matrix file')
    parser.add_argument('fasta_dir', help='Directory containing consensus FASTA files')
    parser.add_argument('output_dir', help='Output directory for prepared FASTA files')
    parser.add_argument('--group1', help='Comma-separated list of samples in group 1')
    parser.add_argument('--group2', help='Comma-separated list of samples in group 2')
    return parser.parse_args()

def load_fasta_sequences(fasta_dir):
    """Load all FASTA sequences from the directory"""
    fasta_dict = {}
    for file in glob.glob(f"{fasta_dir}/*.consensus.fa"):
        with open(file) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                fasta_dict[record.id] = str(record.seq)
    return fasta_dict

def classify_ortholog_groups(df, group1, group2):
    """
    Classify ortholog groups into different categories:
    - monomorphic_both_groups: both groups are monomorphic (all present or all absent)
    - polymorphic_group1: group1 is polymorphic (some present, some absent)
    - other: any other combination
    """
    classification = {}
    
    # Get all unique ortholog IDs
    ortholog_ids = df['ortholog_id'].unique()
    
    for ortholog_id in ortholog_ids:
        # Filter rows for this ortholog
        ortholog_df = df[df['ortholog_id'] == ortholog_id]
        
        # Determine presence/absence for each group
        g1_rows = ortholog_df[ortholog_df['sample'].isin(group1)]
        g2_rows = ortholog_df[ortholog_df['sample'].isin(group2)]
        
        # Get presence status (1=present, 2=absent, 3=missing)
        g1_presence = g1_rows['presence'].to_list()
        g2_presence = g2_rows['presence'].to_list()
        
        # Skip if no data for either group
        if not g1_presence or not g2_presence:
            continue
            
        # Skip if all missing data
        if all(p == 3 for p in g1_presence) and all(p == 3 for p in g2_presence):
            continue
            
        # Skip if any missing data (3)
        if 3 in g1_presence or 3 in g2_presence:
            continue
        
        # Count occurrences of present (1) and absent (2) in each group
        g1_present_count = g1_presence.count(1)
        g1_absent_count = g1_presence.count(2)
        g2_present_count = g2_presence.count(1)
        g2_absent_count = g2_presence.count(2)
        
        # Check if groups are monomorphic (all 1's or all 2's)
        g1_monomorphic = g1_present_count == len(g1_presence) or g1_absent_count == len(g1_presence)
        g2_monomorphic = g2_present_count == len(g2_presence) or g2_absent_count == len(g2_presence)
        
        # Classify the ortholog group
        if g1_monomorphic and g2_monomorphic:
            category = "monomorphic_both_groups"
            g1_status = "present" if g1_present_count == len(g1_presence) else "absent"
            g2_status = "present" if g2_present_count == len(g2_presence) else "absent"
            
            # Include groups even if both are absent (for diversity analysis)
            # Previously skipped groups where both were absent, but we want to analyze those too
                
        elif not g1_monomorphic and g2_monomorphic:
            category = "polymorphic_group1"
            g2_status = "present" if g2_present_count == len(g2_presence) else "absent"
            
            # Skip if there's only one present or one absent sample in group1
            if g1_present_count <= 1 or g1_absent_count <= 1:
                continue
        else:
            # Skip other cases
            continue
        
        # Store classification information
        classification[ortholog_id] = {
            "category": category,
            "g1_status": g1_status if category == "monomorphic_both_groups" else None,
            "g2_status": g2_status,
            "g1_present_count": g1_present_count,
            "g1_absent_count": g1_absent_count,
            "g2_present_count": g2_present_count,
            "g2_absent_count": g2_absent_count
        }
    
    return classification

def prepare_fasta_files(df, fasta_dict, classification, output_dir, group1, group2):
    """Prepare FASTA files for each ortholog group based on its classification"""
    
    for ortholog_id, info in classification.items():
        category = info["category"]
        
        # Filter rows for this ortholog
        ortholog_df = df[df['ortholog_id'] == ortholog_id]
        
        # Process for both left and right flanking regions
        for side in ["left", "right"]:
            if category == "monomorphic_both_groups":
                # Create FASTA for group1
                g1_file = os.path.join(output_dir, f"{ortholog_id}.g1.{side}_flank.fa")
                with open(g1_file, 'w') as f:
                    for _, row in ortholog_df[ortholog_df['sample'].isin(group1)].iterrows():
                        if row['presence'] in [1, 2]:  # Include both present (1) and absent (2) sequences
                            flank_id = get_flank_id(row, side)
                            if flank_id in fasta_dict:
                                f.write(f">{row['sample']}_{flank_id}\n{fasta_dict[flank_id]}\n")
                
                # Create FASTA for group2
                g2_file = os.path.join(output_dir, f"{ortholog_id}.g2.{side}_flank.fa")
                with open(g2_file, 'w') as f:
                    for _, row in ortholog_df[ortholog_df['sample'].isin(group2)].iterrows():
                        if row['presence'] in [1, 2]:  # Include both present (1) and absent (2) sequences
                            flank_id = get_flank_id(row, side)
                            if flank_id in fasta_dict:
                                f.write(f">{row['sample']}_{flank_id}\n{fasta_dict[flank_id]}\n")
            
            elif category == "polymorphic_group1":
                # Create FASTA for present sequences in group1
                present_file = os.path.join(output_dir, f"{ortholog_id}.g1_present.{side}_flank.fa")
                with open(present_file, 'w') as f:
                    for _, row in ortholog_df[(ortholog_df['sample'].isin(group1)) & (ortholog_df['presence'] == 1)].iterrows():
                        flank_id = get_flank_id(row, side)
                        if flank_id in fasta_dict:
                            f.write(f">{row['sample']}_{flank_id}\n{fasta_dict[flank_id]}\n")
                
                # Create FASTA for absent sequences in group1 (if any have sequences)
                absent_file = os.path.join(output_dir, f"{ortholog_id}.g1_absent.{side}_flank.fa")
                with open(absent_file, 'w') as f:
                    for _, row in ortholog_df[(ortholog_df['sample'].isin(group1)) & (ortholog_df['presence'] == 2)].iterrows():
                        # For absent sequences, we need to check if they have flanking sequences
                        if row['start'] != -1 and row['end'] != -1:
                            flank_id = get_flank_id(row, side)
                            if flank_id in fasta_dict:
                                f.write(f">{row['sample']}_{flank_id}\n{fasta_dict[flank_id]}\n")

def get_flank_id(row, side):
    """Get the flank ID for a given row and side"""
    if row['start'] == -1 or row['end'] == -1:
        return None
        
    start = int(row['start'])
    end = int(row['end'])
    contig = row['contig']
    orientation = row['orientation']
    
    # Handle reverse orientation by swapping biological left/right flanks
    if orientation == 'reverse':
        if side == "left":  # Biological 5' flank
            flank_start, flank_end = end - 100, end
        else:  # Biological 3' flank  
            flank_start, flank_end = start, start + 100
    else:  # Forward orientation
        if side == "left":
            flank_start, flank_end = start, start + 100
        else:
            flank_start, flank_end = end - 100, end
        
    return f"{contig}:{flank_start}-{flank_end}"

def main():
    args = parse_arguments()
    
    # Parse group lists
    group1 = args.group1.split(',')
    group2 = args.group2.split(',')
    
    print(f"Reading genotype matrix from {args.genotype_matrix}")
    df = pd.read_csv(args.genotype_matrix, sep="\t")
    
    print(f"Loading FASTA sequences from {args.fasta_dir}")
    fasta_dict = load_fasta_sequences(args.fasta_dir)
    
    print("Classifying ortholog groups")
    classification = classify_ortholog_groups(df, group1, group2)
    
    print(f"Preparing FASTA files for {len(classification)} ortholog groups")
    prepare_fasta_files(df, fasta_dict, classification, args.output_dir, group1, group2)
    
    # Save classification dictionary
    classification_file = os.path.join(args.output_dir, "ortholog_classification.json")
    with open(classification_file, "w") as f:
        json.dump(classification, f, indent=2)
    
    print(f"Classification saved to {classification_file}")
    print("Done!")

if __name__ == "__main__":
    main() 