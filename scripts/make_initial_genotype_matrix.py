import pandas as pd
import os
from collections import defaultdict
import time 
import sys

pd.options.mode.chained_assignment = None

def process_sample_file(tsv: str) -> pd.DataFrame:
    
    # Load the original dataframe
    df = pd.read_csv(tsv, sep="\t")
            
    return df.sort_values(['ortholog_id', 'sample'], 
                              ignore_index=True)

def create_location_index(df: pd.DataFrame) -> dict:
    """Create an index mapping contig+start to ortholog_id"""
    location_index = defaultdict(str)
    # Only index non-missing rows
    for _, row in df[df['presence'] != 3].iterrows():
        key = (row['contig'], row['left_flank_start'])
        location_index[key] = row['ortholog_id']
    return location_index

def count_missing(group: pd.DataFrame) -> int:
    """Count number of missing (3) values in a group"""
    return (group['presence'] == 3).sum()

def combine_ortholog_files(sample_dfs: dict) -> pd.DataFrame:
    """Combine ortholog groups from multiple samples"""
    if not sample_dfs:
        return pd.DataFrame()
    
    # Start with first sample's df as final df
    final_df = sample_dfs['CCMP1545'].copy()
    sample_dfs = {sample : sample_dfs[sample] for sample in sample_dfs if sample != 'CCMP1545'}
    # Process each subsequent sample's df
    for sample in sample_dfs.keys():
        current_df = sample_dfs[sample]
        timestamp = time.strftime("%Y-%m-%d %H:%M:%S")  # Format: YYYY-MM-DD HH:MM:SS
        print(f"next sample {sample}\t{timestamp}")
        new_final_df = final_df.copy()
        
        # Create index for current_df
        current_index = create_location_index(current_df)
        # Track which groups in current_df were matched
        matched_groups = set()
        
        # Check each ortholog group in final_df
        for group_id in final_df['ortholog_id'].unique():
            final_group = final_df[final_df['ortholog_id'] == group_id]
            found_match = False
            
            # Check each non-missing row in the group
            for _, row in final_group[final_group['presence'] != 3].iterrows():
                key = (row['contig'], row['left_flank_start'])
                
                # If we find a match in the index
                if key in current_index:
                    # print("group matched", key, current_index[key])
                    matching_group_id = current_index[key]
                    matched_groups.add(matching_group_id)
                    curr_group = current_df[current_df['ortholog_id'] == matching_group_id]
                    
                    # Compare missing data
                    if count_missing(curr_group) < count_missing(final_group):
                        # Remove old group from new_final_df
                        new_final_df = new_final_df[new_final_df['ortholog_id'] != group_id]
                        # Add new group with original group_id
                        replacement_group = curr_group.copy()
                        replacement_group['ortholog_id'] = group_id
                        new_final_df = pd.concat([new_final_df, replacement_group])
                        found_match = True
                        break
                
                if found_match:
                    break
                    
        # Find unmatched groups
        unmatched_groups = set(current_df['ortholog_id'].unique()) - matched_groups
        if unmatched_groups:
            print(f"Found {len(unmatched_groups)} unmatched groups in {sample}")
            # Generate new ortholog IDs for unmatched groups
            max_id = int(new_final_df['ortholog_id'].str.extract(r'ortholog_id_(\d+)')[0].astype(int).max())
            
            for i, group_id in enumerate(unmatched_groups, 1):
                # Get unmatched group rows
                unmatched_group = current_df[current_df['ortholog_id'] == group_id]
                # Assign new ortholog ID
                new_id = f'ortholog_id_{max_id + i:04d}'
                unmatched_group['ortholog_id'] = new_id
                # Add to new_final_df
                new_final_df = pd.concat([new_final_df, unmatched_group])
        
        final_df = new_final_df
    
    return final_df

def main():
    # Define all samples
    samples = ["CCMP1545", "RCC1749", "CCMP490", "RCC114", "RCC1614", "RCC1698", 
               "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", 
               "RCC692", "RCC693", "RCC833", "RCC835"]
    # samples = ["CCMP1545", "RCC1749", "CCMP490", "RCC3052", "RCC1614"]
    # Set paths
    input_dir = "/scratch1/chris/introner-genotyping-pipeline/graph_blast_results"
    output_file = sys.argv[1]
    
    # Process each sample file
    sample_dfs = {}
    for sample in samples:
        tsv_path = os.path.join(input_dir, f"{sample}.ortholog_results.adjusted.tsv")
        if os.path.exists(tsv_path):
            sample_dfs[sample] = process_sample_file(tsv_path)
    
    # Combine ortholog groups
    final_df = combine_ortholog_files(sample_dfs)
    
    # Save results
    final_df.sort_values(['ortholog_id', 'sample']).to_csv(
        output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")  # Format: YYYY-MM-DD HH:MM:SS
    print(f"Done {timestamp}")