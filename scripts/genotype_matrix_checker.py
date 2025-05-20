import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def analyze_genotype_matrix(file_path):
    """Compute and display statistics for an introner genotype matrix."""
    # Read the matrix
    df = pd.read_csv(file_path, sep='\t')
    print(f"Loaded matrix with {len(df)} rows")
    
    # Basic counts
    ortholog_count = df['ortholog_id'].nunique()
    sample_count = df['sample'].nunique()
    print(f"\nTotal ortholog groups: {ortholog_count}")
    print(f"Total samples: {sample_count}")
    
    # Convert to wide format for easier analysis
    print("\nConverting to wide format for analysis...")
    pivot_df = df.pivot(index='ortholog_id', columns='sample', values='presence')
    
    # Count presence values
    print("\n--- Presence Value Distribution ---")
    presence_counts = df['presence'].value_counts().sort_index()
    total = len(df)
    for value, count in presence_counts.items():
        print(f"Presence {value}: {count} ({count/total*100:.1f}%)")
    
    # Ortholog group statistics
    print("\n--- Ortholog Group Statistics ---")
    
    # All samples have introner present (all 1's)
    all_present = (pivot_df == 1).all(axis=1).sum()
    print(f"Ortholog groups where all samples have introner present: {all_present} ({all_present/ortholog_count*100:.1f}%)")
    
    # All samples have introner absent (all 2's)
    all_absent = (pivot_df == 2).all(axis=1).sum()
    print(f"Ortholog groups where all samples have introner absent: {all_absent} ({all_absent/ortholog_count*100:.1f}%)")
    
    # Some present, some absent (mix of 1's and 2's, no 3's)
    polymorphic = ((pivot_df == 1).any(axis=1) & 
                  (pivot_df == 2).any(axis=1) & 
                  ~(pivot_df == 3).any(axis=1)).sum()
    print(f"Polymorphic ortholog groups (mix of present/absent, no missing data): {polymorphic} ({polymorphic/ortholog_count*100:.1f}%)")
    
    # Groups with missing data
    has_missing = (pivot_df == 3).any(axis=1).sum()
    print(f"Ortholog groups with at least one missing value: {has_missing} ({has_missing/ortholog_count*100:.1f}%)")
    
    # All missing (all 3's)
    all_missing = (pivot_df == 3).all(axis=1).sum()
    print(f"Ortholog groups where all samples have missing data: {all_missing} ({all_missing/ortholog_count*100:.1f}%)")
    
    # Per-sample statistics
    print("\n--- Per-Sample Statistics ---")
    sample_stats = pd.DataFrame({
        'Present (1)': (pivot_df == 1).sum(),
        'Absent (2)': (pivot_df == 2).sum(),
        'Missing (3)': (pivot_df == 3).sum()
    }).reset_index()
    sample_stats['Total'] = sample_stats['Present (1)'] + sample_stats['Absent (2)'] + sample_stats['Missing (3)']
    sample_stats['Present %'] = sample_stats['Present (1)'] / sample_stats['Total'] * 100
    sample_stats['Absent %'] = sample_stats['Absent (2)'] / sample_stats['Total'] * 100
    sample_stats['Missing %'] = sample_stats['Missing (3)'] / sample_stats['Total'] * 100
    
    print(sample_stats)
    
    # Generate additional useful statistics
    print("\n--- Additional Statistics ---")
    
    # Family distribution
    if 'family' in df.columns and not df['family'].isnull().all():
        print("\nTop 10 most common families:")
        family_counts = df[df['family'] != '']['family'].value_counts().head(10)
        for family, count in family_counts.items():
            print(f"  {family}: {count}")
    
    # Gene distribution
    if 'gene' in df.columns and not df['gene'].isnull().all():
        print("\nTop 10 most common genes:")
        gene_counts = df[df['gene'] != '']['gene'].value_counts().head(10)
        for gene, count in gene_counts.items():
            print(f"  {gene}: {count}")
    
    # Calculate allele frequency spectrum for polymorphic groups
    polymorphic_matrix = pivot_df[((pivot_df == 1).any(axis=1) & (pivot_df == 2).any(axis=1))]
    if not polymorphic_matrix.empty:
        print("\nAllele Frequency Spectrum (how many groups have X samples with presence):")
        afs = (polymorphic_matrix == 1).sum(axis=1).value_counts().sort_index()
        for presence_count, num_groups in afs.items():
            print(f"  {presence_count} sample(s) with introner present: {num_groups} ortholog groups")
    
    # Check for duplicate sequences
    check_duplicate_sequences(df)
    
    return df, pivot_df

def check_duplicate_sequences(df):
    """
    Check for duplicate sequences in the genotype matrix.
    A sequence is considered a duplicate if it has the same sample, contig, and coordinates.
    """
    print("\n--- Duplicate Sequence Check ---")
    
    # Filter out rows with missing data
    valid_df = df[df['presence'] != 3].copy()
    
    # Create a unique sequence identifier
    valid_df['seq_id'] = valid_df.apply(
        lambda row: f"{row['sample']}|{row['contig']}:{row['start']}-{row['end']}" 
        if pd.notnull(row['contig']) and pd.notnull(row['start']) and pd.notnull(row['end']) 
        else None, 
        axis=1
    )
    
    # Filter out rows without a valid seq_id
    valid_df = valid_df[valid_df['seq_id'].notnull()]
    
    # Check for duplicates across ortholog groups
    seq_to_ortholog = defaultdict(list)
    for _, row in valid_df.iterrows():
        seq_to_ortholog[row['seq_id']].append(row['ortholog_id'])
    
    # Find sequences that appear in multiple ortholog groups
    duplicates = {seq_id: orthologs for seq_id, orthologs in seq_to_ortholog.items() if len(orthologs) > 1}
    
    if duplicates:
        print(f"Found {len(duplicates)} sequences that appear in multiple ortholog groups:")
        for seq_id, orthologs in duplicates.items():
            parts = seq_id.split('|')
            sample = parts[0]
            coords = parts[1] if len(parts) > 1 else "unknown"
            print(f"  Sequence {sample} at {coords} appears in {len(orthologs)} ortholog groups: {', '.join(orthologs)}")
            
            # Get additional metadata for this sequence
            for _, row in valid_df[valid_df['seq_id'] == seq_id].iterrows():
                print(f"    In {row['ortholog_id']}: presence={row['presence']}, gene={row.get('gene', 'N/A')}, family={row.get('family', 'N/A')}")
            print()
    else:
        print("No duplicate sequences found across ortholog groups!")
    
    # Check for duplicate samples within ortholog groups
    ortho_to_samples = defaultdict(list)
    for _, row in valid_df.iterrows():
        ortho_to_samples[row['ortholog_id']].append(row['sample'])
    
    # Find ortholog groups with duplicate samples
    duplicate_samples = {oid: samples for oid, samples in ortho_to_samples.items() 
                         if len(samples) != len(set(samples))}
    
    if duplicate_samples:
        print(f"\nFound {len(duplicate_samples)} ortholog groups with duplicate samples:")
        for oid, samples in duplicate_samples.items():
            # Count sample occurrences
            sample_counts = {}
            for sample in samples:
                sample_counts[sample] = sample_counts.get(sample, 0) + 1
            
            # Report duplicates
            duplicates = {s: c for s, c in sample_counts.items() if c > 1}
            print(f"  Ortholog group {oid} has duplicate samples: {duplicates}")
    else:
        print("\nNo duplicate samples found within ortholog groups!")

    # Check for proximity-based potential duplicates
    print("\nChecking for proximity-based potential duplicates...")
    proximity_duplicates = check_proximity_duplicates(valid_df, window=100)
    
    if proximity_duplicates:
        print(f"Found {len(proximity_duplicates)} potential proximity-based duplicates:")
        for (sample, contig), seq_groups in proximity_duplicates.items():
            print(f"  Sample {sample}, contig {contig} has {len(seq_groups)} sequence clusters that may be duplicates:")
            for i, group in enumerate(seq_groups, 1):
                print(f"    Cluster {i}: {len(group)} sequences")
                for seq in group[:5]:  # Show at most 5 sequences per cluster
                    # Get ortholog IDs for this sequence
                    orthologs = [row['ortholog_id'] for _, row in valid_df[
                        (valid_df['sample'] == sample) & 
                        (valid_df['contig'] == contig) & 
                        (valid_df['start'] == seq[0]) & 
                        (valid_df['end'] == seq[1])
                    ].iterrows()]
                    print(f"      {seq[0]}-{seq[1]} in ortholog groups: {', '.join(orthologs)}")
                if len(group) > 5:
                    print(f"      ... and {len(group) - 5} more")
    else:
        print("No proximity-based duplicates found!")

def check_proximity_duplicates(df, window=100):
    """
    Check for sequences from the same sample and contig that are very close to each other,
    which may indicate potential duplicates.
    
    Args:
        df: DataFrame with sequence information
        window: Maximum distance between sequences to consider them potential duplicates
        
    Returns:
        Dictionary of potential duplicates
    """
    # Group by sample and contig
    proximity_duplicates = {}
    
    # Get unique sample-contig combinations
    sample_contigs = df[['sample', 'contig']].drop_duplicates()
    
    for _, row in sample_contigs.iterrows():
        sample = row['sample']
        contig = row['contig']
        
        # Get all sequences for this sample and contig
        sample_seqs = df[(df['sample'] == sample) & (df['contig'] == contig)].copy()
        
        if len(sample_seqs) <= 1:
            continue
        
        # Sort by start position
        sample_seqs = sample_seqs.sort_values('start')
        
        # Extract start and end positions
        positions = [(int(row['start']), int(row['end'])) 
                     for _, row in sample_seqs.iterrows() 
                     if pd.notnull(row['start']) and pd.notnull(row['end'])]
        
        if not positions:
            continue
            
        # Find clusters of sequences that are close to each other
        clusters = []
        current_cluster = [positions[0]]
        
        for i in range(1, len(positions)):
            prev_end = positions[i-1][1]
            curr_start = positions[i][0]
            
            if curr_start - prev_end <= window:
                # This sequence is close to the previous one, add to current cluster
                current_cluster.append(positions[i])
            else:
                # Start a new cluster
                if len(current_cluster) > 1:
                    clusters.append(current_cluster)
                current_cluster = [positions[i]]
        
        # Add the last cluster if it has multiple sequences
        if len(current_cluster) > 1:
            clusters.append(current_cluster)
        
        # Store clusters if any were found
        if clusters:
            proximity_duplicates[(sample, contig)] = clusters
    
    return proximity_duplicates

def plot_statistics(df, pivot_df, output_prefix="genotype_stats"):
    """Generate plots for the genotype matrix statistics."""
    # Create a directory for plots if it doesn't exist
    import os
    if not os.path.exists("plots"):
        os.makedirs("plots")
    
    # Presence distribution
    plt.figure(figsize=(10, 6))
    presence_counts = df['presence'].value_counts().sort_index()
    sns.barplot(x=presence_counts.index, y=presence_counts.values)
    plt.xlabel('Presence Value')
    plt.ylabel('Count')
    plt.title('Distribution of Presence Values')
    plt.xticks([0, 1, 2], ['1 (Present)', '2 (Absent)', '3 (Missing)'])
    plt.tight_layout()
    plt.savefig(f"plots/{output_prefix}_presence_distribution.png")
    
    # Heatmap of presence/absence
    if len(pivot_df) <= 100:  # Only generate heatmap for small matrices
        plt.figure(figsize=(12, 10))
        sns.heatmap(pivot_df, cmap='viridis', cbar_kws={'label': 'Presence Value'})
        plt.title('Genotype Matrix Heatmap')
        plt.tight_layout()
        plt.savefig(f"plots/{output_prefix}_heatmap.png")
    
    # Allele frequency spectrum
    polymorphic_matrix = pivot_df[((pivot_df == 1).any(axis=1) & (pivot_df == 2).any(axis=1))]
    if not polymorphic_matrix.empty:
        plt.figure(figsize=(10, 6))
        afs = (polymorphic_matrix == 1).sum(axis=1).value_counts().sort_index()
        sns.barplot(x=afs.index, y=afs.values)
        plt.xlabel('Number of Samples with Introner Present')
        plt.ylabel('Number of Ortholog Groups')
        plt.title('Allele Frequency Spectrum')
        plt.tight_layout()
        plt.savefig(f"plots/{output_prefix}_afs.png")
    
    # Per-sample statistics
    sample_presence = (pivot_df == 1).sum()
    sample_absence = (pivot_df == 2).sum()
    sample_missing = (pivot_df == 3).sum()
    
    plt.figure(figsize=(12, 6))
    x = range(len(sample_presence))
    width = 0.3
    plt.bar([i-width for i in x], sample_presence, width=width, label='Present')
    plt.bar(x, sample_absence, width=width, label='Absent')
    plt.bar([i+width for i in x], sample_missing, width=width, label='Missing')
    plt.xlabel('Sample')
    plt.ylabel('Count')
    plt.title('Presence/Absence/Missing by Sample')
    plt.xticks(x, sample_presence.index, rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"plots/{output_prefix}_sample_distribution.png")
    
    print(f"\nPlots saved to plots/{output_prefix}_*.png")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python analyze_genotype_matrix.py <genotype_matrix.tsv> [output_prefix]")
        sys.exit(1)
    
    file_path = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else "genotype_stats"
    
    df, pivot_df = analyze_genotype_matrix(file_path)
    plot_statistics(df, pivot_df, output_prefix)