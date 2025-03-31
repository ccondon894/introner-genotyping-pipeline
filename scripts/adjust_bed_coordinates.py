#!/usr/bin/env python3

import sys
import os
import pandas as pd

def process_df(df, sample):
    """Process a single sample's ortholog_results.tsv file with vectorized operations"""
    samples = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", 
           "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", 
           "RCC692", "RCC693", "RCC833", "RCC835"]

    # Extract contig, start, and end directly using vectorized string operations
    new_df = df[['ortholog_id']].copy()
    new_df = new_df.drop_duplicates(ignore_index=True)
    new_df['contig'] = new_df['ortholog_id'].str.split(":").str[0]
    coords = new_df['ortholog_id'].str.split(":").str[1].str.split("|").str[0]
    new_df['left_flank_start'] = coords.str.split("-").str[0].astype(int)
    new_df['right_flank_end'] = coords.str.split("-").str[1].astype(int)
    new_df['presence'] = 0  # Assuming presence is always 0
    new_df['sample'] = sample

    # Define column order explicitly to avoid any misalignment
    new_df = new_df[['ortholog_id', 'contig', 'left_flank_start', 'right_flank_end', 'presence', 'sample']]

    # Combine the old and new dataframes efficiently
    combined_df = pd.concat([df, new_df], ignore_index=True)
    combined_df = combined_df.sort_values(['ortholog_id', 'sample'])
    # Assign generic ortholog_id values
    unique_ids = combined_df['ortholog_id'].unique()
    id_mapping = {orth_id: f'ortholog_id_{i:04d}' for i, orth_id in enumerate(unique_ids, start=1)}
    combined_df['ortholog_id'] = combined_df['ortholog_id'].map(id_mapping)

    # List of expected samples
    samples = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", 
            "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", 
            "RCC692", "RCC693", "RCC833", "RCC835"]

    # Identify unique ortholog IDs
    unique_orthologs = combined_df['ortholog_id'].unique()

    # Group by ortholog_id to find missing samples
    existing_samples = combined_df.groupby('ortholog_id')['sample'].apply(set)

    # List to store missing sample rows
    missing_rows = []

    for ortholog_id in unique_orthologs:
        # Find missing samples for this ortholog_id
        present_samples = existing_samples.get(ortholog_id, set())
        missing_samples = set(samples) - present_samples
        
        # Create rows for missing samples
        for sample in missing_samples:
            missing_rows.append({
                'ortholog_id': ortholog_id,
                'contig': '?',
                'left_flank_start': 0,
                'right_flank_end': 0,
                'presence': 3,
                'sample': sample
            })

    # Convert missing rows into a dataframe
    missing_df = pd.DataFrame(missing_rows)

    # Combine with the original dataframe
    final_df = pd.concat([combined_df, missing_df], ignore_index=True)
    # Update rows where presence is 3 to have coordinates of 0
    final_df.loc[final_df['presence'] == 3, ['left_flank_start', 'right_flank_end']] = 0
    final_df.loc[final_df['presence'] == 3, ['contig']] = '?'
    final_df = final_df.astype({
        'left_flank_start': 'Int64',
        'right_flank_end': 'Int64',
        'presence': 'Int64'})
    
    return final_df.sort_values(['ortholog_id', 'sample'], 
                              ignore_index=True)

def overlap_size(start1, end1, start2, end2):
    """
    Returns the number of overlapping bases between two intervals.
    Overlap = (min(end1, end2) - max(start1, start2) + 1)
    Returns 0 if there is no overlap.
    """
    return max(0, min(end1, end2) - max(start1, start2) + 1)

def find_best_bed_match(row, df_bed):
    """
    Given a single 'results' row and a DataFrame of BED intervals (for the same sample),
    find the BED interval with the largest overlap. Returns that BED interval row or
    None if there is no overlapping interval with at least 90% overlap.
    """
    contig = row['contig']
    r_start = row['left_flank_start']
    r_end = row['right_flank_end']
    
    # Check for NaN values in input
    if pd.isna(r_start) or pd.isna(r_end) or pd.isna(contig):
        return None
        
    # Check for valid coordinates
    query_length = r_end - r_start
    if query_length <= 0:
        return None

    # Filter to the same contig in the bed
    df_same_contig = df_bed[df_bed['contig'] == contig]
    if df_same_contig.empty:
        return None

    # Compute overlap for each row in df_same_contig
    overlaps = df_same_contig.apply(
        lambda bed_row: overlap_size(r_start, r_end, bed_row['start'], bed_row['end']),
        axis=1
    )
    
    # Check if there are any non-zero overlaps
    if (overlaps == 0).all():
        return None

    # Calculate overlap percentages
    overlap_percentages = overlaps / query_length * 100
    
    # If there's no overlap > 90%, return None
    if overlap_percentages.max() < 90:
        return None

    # Find the index of the maximum overlap that meets the threshold
    max_overlap_index = overlap_percentages.idxmax()
    return df_same_contig.loc[max_overlap_index]

def adjust_results_coordinates(results_file, bed_folder):
    """
    Adjusts the coordinates in `results_file` by matching them to intervals in
    the corresponding BED files located in `bed_folder`.
    If no matching bed region is found, keeps the original coordinates.

    The results file should have at least columns:
      [contig, left_flank_start, right_flank_end, sample]

    Each BED file is assumed to be named <sample>.bed, with columns [contig, start, end].
    """
    # Load the results file
    df_results = pd.read_csv(results_file, sep="\t")

    # Prepare a list to hold all rows across all samples
    all_rows = []

    # Group results by sample
    grouped = df_results.groupby('sample')
    
    for sample, group_df in grouped:
        # Path to the sample's BED file
        bed_path = os.path.join(bed_folder, f"{sample}.candidate_loci_plus_flanks.bed")
        
        if not os.path.exists(bed_path):
            print(f"[Warning] BED file not found for sample '{sample}' at {bed_path}. Keeping original coordinates.")
            # Keep all rows for this sample as-is
            all_rows.append(group_df)
            continue
        
        # Load the BED file for this sample
        df_bed = pd.read_csv(bed_path, sep=" ", names=['contig','start','end', 'seq_id'], header=None)

        # We'll store all rows, updated or original
        group_rows = []
        
        # Process each row in this group
        for idx, row in group_df.iterrows():
            best_match = find_best_bed_match(row, df_bed)
            if best_match is not None:
                # Create an updated copy of the row
                updated_row = row.copy()
                updated_row['left_flank_start'] = best_match['start']
                updated_row['right_flank_end']  = best_match['end']
                group_rows.append(updated_row)
            else:
                # No overlap found -> keep original row
                group_rows.append(row)
        
        # Convert list of rows to a DataFrame and append
        group_df = pd.DataFrame(group_rows)
        all_rows.append(group_df)
    
    # Concatenate all rows
    df_final = pd.concat(all_rows, ignore_index=True)
    return df_final

def main():
    """
    Usage:
        python adjust_coordinates.py <results_file> <bed_folder> <output_file>

    Example:
        python adjust_coordinates.py my_results.tsv bed_files adjusted_results.tsv
    """
    if len(sys.argv) != 4:
        print("Usage: python adjust_coordinates.py <results_file> <bed_folder> <output_file>")
        sys.exit(1)

    results_file = sys.argv[1]
    bed_folder = sys.argv[2]
    output_file = sys.argv[3]

    sample = results_file.split("/")[-1].split(".")[0]
    # Adjust coordinates
    df_adjusted = adjust_results_coordinates(results_file, bed_folder)

    
    # Then process the adjusted results with process_sample_file
    df_final = process_df(df_adjusted, sample)
    
    # Save to output file
    df_final.to_csv(output_file, sep="\t", index=False)
    
    
    print(f"Done adjusting results. Saved to: {output_file}")

if __name__ == "__main__":
    main()