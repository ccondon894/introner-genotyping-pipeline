#!/usr/bin/env python3

import argparse
import os
import glob
import json
import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description='Calculate diversity metrics (Pi and dXY) from alignments')
    parser.add_argument('--alignment_dir', help='Directory containing alignment files')
    parser.add_argument('--classification', help='Path to ortholog classification JSON file')
    parser.add_argument('--output', help='Output file for diversity metrics')
    parser.add_argument('--verbose', action='store_true', help='Print verbose output')
    return parser.parse_args()

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

def calculate_between_groups_dxy(sequences1, sequences2):
    """Calculate average dXY between two groups of sequences"""
    if not sequences1 or not sequences2:
        return 0
    
    total_dxy = 0
    comparisons = 0
    
    for seq1 in sequences1:
        for seq2 in sequences2:
            dxy = calculate_pairwise_pi(seq1, seq2)
            total_dxy += dxy
            comparisons += 1
    
    return total_dxy / comparisons if comparisons > 0 else 0

def process_monomorphic_groups(alignment_dir, classification, verbose=False):
    """Process ortholog groups where both groups are monomorphic"""
    # Store results by ortholog_id
    results_by_ortholog = {}
    
    for ortholog_id, info in classification.items():
        if info["category"] != "monomorphic_both_groups":
            continue
            
        g1_status = info["g1_status"]
        g2_status = info["g2_status"]
        
        # Initialize side-specific results
        side_results = {}
        
        for side in ["left", "right"]:
            # Get alignment files
            g1_file = os.path.join(alignment_dir, f"{ortholog_id}.g1.{side}_flank.mafft.fa")
            g2_file = os.path.join(alignment_dir, f"{ortholog_id}.g2.{side}_flank.mafft.fa")
            combined_file = os.path.join(alignment_dir, f"{ortholog_id}.g1_g2_combined.{side}_flank.mafft.fa")
            
            # Check which files need to exist based on presence/absence state
            files_to_check = []
            if g1_status == "present" or g1_status == "absent":
                files_to_check.append(g1_file)
            if g2_status == "present" or g2_status == "absent":
                files_to_check.append(g2_file)
            
            # Only check combined file if both groups have sequences
            if (g1_status == "present" or g1_status == "absent") and (g2_status == "present" or g2_status == "absent"):
                files_to_check.append(combined_file)
            
            # Skip if required files don't exist
            missing_files = [f for f in files_to_check if not os.path.exists(f)]
            if missing_files:
                if verbose:
                    print(f"Skipping {ortholog_id} {side}: Missing files: {', '.join(missing_files)}")
                continue
                
            try:
                # Initialize metrics
                g1_pi = None
                g2_pi = None
                dxy = None
                
                # Calculate metrics based on presence/absence state
                
                # Process Group 1 if it has sequences
                if g1_status == "present" or g1_status == "absent":
                    g1_seqs = [str(record.seq) for record in SeqIO.parse(g1_file, "fasta")]
                    if len(g1_seqs) > 1:
                        g1_pi = calculate_group_pi(g1_seqs)
                
                # Process Group 2 if it has sequences
                if g2_status == "present" or g2_status == "absent":
                    g2_seqs = [str(record.seq) for record in SeqIO.parse(g2_file, "fasta")]
                    if len(g2_seqs) > 1:
                        g2_pi = calculate_group_pi(g2_seqs)
                
                # Calculate dXY only if both groups have sequences
                if (g1_status == "present" or g1_status == "absent") and (g2_status == "present" or g2_status == "absent"):
                    # Load combined sequences for dXY calculation
                    combined_records = list(SeqIO.parse(combined_file, "fasta"))
                    
                    # Get g1 sequence IDs to separate them in the combined alignment
                    g1_ids = set(record.id for record in SeqIO.parse(g1_file, "fasta"))
                    
                    # Separate sequences from the combined alignment
                    g1_combined_seqs = [str(record.seq) for record in combined_records if record.id in g1_ids]
                    g2_combined_seqs = [str(record.seq) for record in combined_records if record.id not in g1_ids]
                    
                    if g1_combined_seqs and g2_combined_seqs:
                        # Calculate dXY between groups using sequences from the combined alignment
                        dxy = calculate_between_groups_dxy(g1_combined_seqs, g2_combined_seqs)
                
                # Store results for this side
                side_results[side] = {
                    "g1_pi": g1_pi,
                    "g2_pi": g2_pi,
                    "dxy": dxy
                }
                
            except Exception as e:
                if verbose:
                    print(f"Error processing {ortholog_id} {side}: {e}")
                continue
        
        # Only include ortholog if we have results for at least one side
        if side_results:
            # Calculate average metrics across sides
            g1_pi_values = [data["g1_pi"] for side, data in side_results.items() if data["g1_pi"] is not None]
            g2_pi_values = [data["g2_pi"] for side, data in side_results.items() if data["g2_pi"] is not None]
            dxy_values = [data["dxy"] for side, data in side_results.items() if data["dxy"] is not None]
            
            results_by_ortholog[ortholog_id] = {
                "ortholog_id": ortholog_id,
                "category": "monomorphic_both_groups",
                "g1_status": g1_status,
                "g2_status": g2_status,
                "g1_pi": np.mean(g1_pi_values) if g1_pi_values else None,
                "g2_pi": np.mean(g2_pi_values) if g2_pi_values else None,
                "dxy": np.mean(dxy_values) if dxy_values else None,
                "g1_present_pi": None,
                "g1_absent_pi": None,
                "g1_present_absent_pi": None,
                "sides_processed": ",".join(side_results.keys()),
                "g1_present_count": info["g1_present_count"],
                "g1_absent_count": info["g1_absent_count"],
                "g2_present_count": info["g2_present_count"],
                "g2_absent_count": info["g2_absent_count"]
            }
    
    return list(results_by_ortholog.values())

def process_polymorphic_group1(alignment_dir, classification, verbose=False):
    """Process ortholog groups where group1 is polymorphic"""
    # Store results by ortholog_id
    results_by_ortholog = {}
    
    for ortholog_id, info in classification.items():
        if info["category"] != "polymorphic_group1":
            continue
            
        g2_status = info["g2_status"]
        
        # Check if we have enough sequences in each category
        if info["g1_present_count"] <= 1 or info["g1_absent_count"] <= 1:
            if verbose:
                print(f"Skipping {ortholog_id}: Not enough sequences in present/absent categories")
            continue
        
        # Initialize side-specific results
        side_results = {}
        
        for side in ["left", "right"]:
            # Get alignment files
            present_file = os.path.join(alignment_dir, f"{ortholog_id}.g1_present.{side}_flank.mafft.fa")
            absent_file = os.path.join(alignment_dir, f"{ortholog_id}.g1_absent.{side}_flank.mafft.fa")
            combined_file = os.path.join(alignment_dir, f"{ortholog_id}.g1_combined.{side}_flank.mafft.fa")
            
            # Skip if any file doesn't exist
            if not os.path.exists(present_file):
                if verbose:
                    print(f"Skipping {ortholog_id} {side}: Missing present file: {present_file}")
                continue
                
            if not os.path.exists(absent_file):
                if verbose:
                    print(f"Skipping {ortholog_id} {side}: Missing absent file: {absent_file}")
                continue
                
            if not os.path.exists(combined_file):
                if verbose:
                    print(f"Skipping {ortholog_id} {side}: Missing combined file: {combined_file}")
                continue
                
            try:
                # Load sequences for Pi within group calculations
                present_seqs = [str(record.seq) for record in SeqIO.parse(present_file, "fasta")]
                absent_seqs = [str(record.seq) for record in SeqIO.parse(absent_file, "fasta")]
                
                if not present_seqs:
                    if verbose:
                        print(f"Skipping {ortholog_id} {side}: No present sequences found")
                    continue
                    
                if not absent_seqs:
                    if verbose:
                        print(f"Skipping {ortholog_id} {side}: No absent sequences found")
                    continue
                
                # Calculate Pi within present and absent groups
                present_pi = calculate_group_pi(present_seqs)
                absent_pi = calculate_group_pi(absent_seqs)
                
                # Load combined sequences for Pi between present and absent
                combined_records = list(SeqIO.parse(combined_file, "fasta"))
                
                # Get present sequence IDs to separate them in the combined alignment
                present_ids = set(record.id for record in SeqIO.parse(present_file, "fasta"))
                
                # Separate sequences from the combined alignment
                present_combined_seqs = [str(record.seq) for record in combined_records if record.id in present_ids]
                absent_combined_seqs = [str(record.seq) for record in combined_records if record.id not in present_ids]
                
                if not present_combined_seqs:
                    if verbose:
                        print(f"Skipping {ortholog_id} {side}: No present sequences in combined alignment")
                    continue
                    
                if not absent_combined_seqs:
                    if verbose:
                        print(f"Skipping {ortholog_id} {side}: No absent sequences in combined alignment")
                    continue
                
                # Calculate Pi between present and absent using sequences from the combined alignment
                present_absent_pi = calculate_between_groups_dxy(present_combined_seqs, absent_combined_seqs)
                
                # Store results for this side
                side_results[side] = {
                    "present_pi": present_pi,
                    "absent_pi": absent_pi,
                    "present_absent_pi": present_absent_pi
                }
                
            except Exception as e:
                if verbose:
                    print(f"Error processing {ortholog_id} {side}: {e}")
                continue
        
        # Only include ortholog if we have results for at least one side
        if side_results:
            # Calculate average metrics across sides
            present_pi_values = [data["present_pi"] for side, data in side_results.items()]
            absent_pi_values = [data["absent_pi"] for side, data in side_results.items()]
            present_absent_pi_values = [data["present_absent_pi"] for side, data in side_results.items()]
            
            results_by_ortholog[ortholog_id] = {
                "ortholog_id": ortholog_id,
                "category": "polymorphic_group1",
                "g1_status": "polymorphic",
                "g2_status": g2_status,
                "g1_pi": None,
                "g2_pi": None,
                "dxy": None,
                "g1_present_pi": np.mean(present_pi_values) if present_pi_values else None,
                "g1_absent_pi": np.mean(absent_pi_values) if absent_pi_values else None,
                "g1_present_absent_pi": np.mean(present_absent_pi_values) if present_absent_pi_values else None,
                "sides_processed": ",".join(side_results.keys()),
                "g1_present_count": info["g1_present_count"],
                "g1_absent_count": info["g1_absent_count"],
                "g2_present_count": info["g2_present_count"],
                "g2_absent_count": info["g2_absent_count"]
            }
    
    return list(results_by_ortholog.values())

def main():
    args = parse_arguments()
    
    # Load ortholog classification
    print(f"Loading ortholog classification from {args.classification}")
    with open(args.classification, "r") as f:
        classification = json.load(f)
    
    print(f"Processing {len(classification)} classified ortholog groups")
    
    # Process monomorphic groups
    print("Processing monomorphic groups")
    mono_results = process_monomorphic_groups(args.alignment_dir, classification, args.verbose)
    print(f"Processed {len(mono_results)} monomorphic ortholog groups")
    
    # Process polymorphic group1
    print("Processing polymorphic group1")
    poly_results = process_polymorphic_group1(args.alignment_dir, classification, args.verbose)
    print(f"Processed {len(poly_results)} polymorphic ortholog groups")
    
    # Combine results
    results = mono_results + poly_results
    
    # Convert results to DataFrame and save
    df = pd.DataFrame(results)
    df.to_csv(args.output, sep="\t", index=False)
    
    print(f"Saved diversity metrics for {len(results)} ortholog groups to {args.output}")
    print("Done!")

if __name__ == "__main__":
    main()