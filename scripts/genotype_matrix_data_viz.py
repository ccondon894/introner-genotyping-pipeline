import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

def analyze_genotype_matrix(file_path):
    """Compute and display statistics for an introner genotype matrix with population analysis."""
    # Read the matrix
    df = pd.read_csv(file_path, sep='\t')
    print(f"Loaded matrix with {len(df)} rows")
    
    # Define the two populations
    group2 = ['RCC1749', 'RCC3052']
    group1 = [sample for sample in df['sample'].unique() if sample not in group2]
    
    print(f"\nGroup 1 samples: {', '.join(group1)}")
    print(f"Group 2 samples: {', '.join(group2)}")
    
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
    
    # Population-specific statistics
    print("\n--- Population-Specific Statistics ---")
    
    # Extract group-specific data
    group1_df = pivot_df[group1]
    group2_df = pivot_df[group2]
    
    # Count valid rows (no missing data)
    valid_both_groups = ((group1_df != 3).all(axis=1) & (group2_df != 3).all(axis=1))
    valid_rows = pivot_df[valid_both_groups]
    print(f"Ortholog groups with valid data for both populations: {valid_rows.shape[0]}")
    
    # Calculate population-specific patterns
    group1_all_present = (group1_df == 1).all(axis=1) & valid_both_groups
    group1_all_absent = (group1_df == 2).all(axis=1) & valid_both_groups
    group1_polymorphic = ((group1_df == 1).any(axis=1) & (group1_df == 2).any(axis=1)) & valid_both_groups
    
    group2_all_present = (group2_df == 1).all(axis=1) & valid_both_groups
    group2_all_absent = (group2_df == 2).all(axis=1) & valid_both_groups
    group2_polymorphic = ((group2_df == 1).any(axis=1) & (group2_df == 2).any(axis=1)) & valid_both_groups
    
    print("\nGroup 1 patterns:")
    print(f"  All present: {group1_all_present.sum()} ({group1_all_present.sum()/valid_rows.shape[0]*100:.1f}%)")
    print(f"  All absent: {group1_all_absent.sum()} ({group1_all_absent.sum()/valid_rows.shape[0]*100:.1f}%)")
    print(f"  Polymorphic: {group1_polymorphic.sum()} ({group1_polymorphic.sum()/valid_rows.shape[0]*100:.1f}%)")
    
    print("\nGroup 2 patterns:")
    print(f"  All present: {group2_all_present.sum()} ({group2_all_present.sum()/valid_rows.shape[0]*100:.1f}%)")
    print(f"  All absent: {group2_all_absent.sum()} ({group2_all_absent.sum()/valid_rows.shape[0]*100:.1f}%)")
    print(f"  Polymorphic: {group2_polymorphic.sum()} ({group2_polymorphic.sum()/valid_rows.shape[0]*100:.1f}%)")
    
    # Cross-population patterns
    print("\nCross-population patterns (among valid ortholog groups):")
    g1_present_g2_present = (group1_all_present & group2_all_present).sum()
    g1_present_g2_absent = (group1_all_present & group2_all_absent).sum()
    g1_absent_g2_present = (group1_all_absent & group2_all_present).sum()
    g1_absent_g2_absent = (group1_all_absent & group2_all_absent).sum()
    
    print(f"  Group 1 present, Group 2 present: {g1_present_g2_present} ({g1_present_g2_present/valid_rows.shape[0]*100:.1f}%)")
    print(f"  Group 1 present, Group 2 absent: {g1_present_g2_absent} ({g1_present_g2_absent/valid_rows.shape[0]*100:.1f}%)")
    print(f"  Group 1 absent, Group 2 present: {g1_absent_g2_present} ({g1_absent_g2_present/valid_rows.shape[0]*100:.1f}%)")
    print(f"  Group 1 absent, Group 2 absent: {g1_absent_g2_absent} ({g1_absent_g2_absent/valid_rows.shape[0]*100:.1f}%)")
    
    # Calculate 2D AFS
    afs_2d = calculate_2d_afs(pivot_df, group1, group2)
    
    # Per-sample statistics
    print("\n--- Per-Sample Statistics ---")
    sample_stats = pd.DataFrame({
        'Present (1)': (pivot_df == 1).sum(),
        'Absent (2)': (pivot_df == 2).sum(),
        'Missing (3)': (pivot_df == 3).sum()
    })
    sample_stats['Total'] = sample_stats['Present (1)'] + sample_stats['Absent (2)'] + sample_stats['Missing (3)']
    sample_stats['Present %'] = sample_stats['Present (1)'] / sample_stats['Total'] * 100
    sample_stats['Absent %'] = sample_stats['Absent (2)'] / sample_stats['Total'] * 100
    sample_stats['Missing %'] = sample_stats['Missing (3)'] / sample_stats['Total'] * 100
    
    # Add population info - fixed the error here
    sample_stats['Population'] = [('Group 2' if x in group2 else 'Group 1') for x in sample_stats.index]
    
    # Reset index to make the sample names a column
    sample_stats = sample_stats.reset_index().rename(columns={'index': 'Sample'})
    
    # Group by population
    pop_stats = sample_stats.groupby('Population').agg({
        'Present (1)': 'sum',
        'Absent (2)': 'sum',
        'Missing (3)': 'sum',
        'Total': 'sum'
    })
    pop_stats['Present %'] = pop_stats['Present (1)'] / pop_stats['Total'] * 100
    pop_stats['Absent %'] = pop_stats['Absent (2)'] / pop_stats['Total'] * 100
    pop_stats['Missing %'] = pop_stats['Missing (3)'] / pop_stats['Total'] * 100
    
    print("\nAggregated statistics by population:")
    print(pop_stats)
    
    print("\nDetailed statistics by sample:")
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
    
    return df, pivot_df, group1, group2, afs_2d

def calculate_2d_afs(pivot_df, group1, group2):
    """
    Calculate a 2D allele frequency spectrum for two populations.
    
    Args:
        pivot_df: Pivot table with samples as columns, ortholog groups as rows
        group1: List of sample names in group 1
        group2: List of sample names in group 2
        
    Returns:
        2D numpy array representing the AFS
    """
    # Extract population-specific data
    group1_df = pivot_df[group1]
    group2_df = pivot_df[group2]
    
    # Count only groups where data is complete for both populations
    valid_rows = ((group1_df != 3).all(axis=1) & (group2_df != 3).all(axis=1))
    valid_group1 = group1_df[valid_rows]
    valid_group2 = group2_df[valid_rows]
    
    # Count the number of samples in each group
    n_group1 = len(group1)
    n_group2 = len(group2)
    
    # Initialize 2D AFS matrix
    afs_2d = np.zeros((n_group1 + 1, n_group2 + 1), dtype=int)
    
    # For each valid ortholog group
    for i in range(len(valid_group1)):
        # Count samples with introner present (presence=1) in each group
        count_g1 = (valid_group1.iloc[i] == 1).sum()
        count_g2 = (valid_group2.iloc[i] == 1).sum()
        afs_2d[count_g1, count_g2] += 1
    
    print("\n--- 2D Allele Frequency Spectrum ---")
    print(f"Found {valid_rows.sum()} ortholog groups with complete data for both populations")
    print("AFS matrix (rows: count in Group 1, columns: count in Group 2):")
    print(afs_2d)
    
    return afs_2d

def plot_introner_family_distribution(df, output_prefix="genotype_stats", use_proportions=False):
    """
    Generate bar charts showing counts or proportions of present introners by family across samples.
    Families are on the x-axis with samples grouped for each family.
    
    Args:
        df: The genotype matrix dataframe
        output_prefix: Prefix for output files
        use_proportions: If True, show proportion of each family relative to total present introners in each sample
                        If False, show absolute counts (default)
    """
    # Check if family column exists and has data
    if 'family' not in df.columns or df['family'].isnull().all() or (df['family'] == '').all():
        print("No family information available in the dataset. Skipping family distribution plot.")
        return
    
    # Filter out missing data (presence == 3) and rows without family information
    filtered_df = df[(df['presence'] != 3) & (df['family'] != '') & (~df['family'].isnull())]
    
    if filtered_df.empty:
        print("No valid data available for family distribution plot after filtering.")
        return
    
    # Get top N families by frequency for cleaner visualization
    top_n = 10
    top_families = filtered_df['family'].value_counts().head(top_n).index.tolist()
    
    # Filter data to include only top families
    family_df = filtered_df[filtered_df['family'].isin(top_families)]
    
    # Get unique samples
    samples = sorted(family_df['sample'].unique())
    
    # Calculate counts or proportions
    result_data = {}
    
    # First calculate total present introners per sample (for proportions)
    sample_totals = {}
    for sample in samples:
        # Filter for current sample and presence=1
        sample_data = filtered_df[(filtered_df['sample'] == sample) & (filtered_df['presence'] == 1)]
        # Store total count of present introners for this sample
        sample_totals[sample] = len(sample_data)
    
    # Now calculate values for each family
    for family in top_families:
        family_data = family_df[family_df['family'] == family]
        sample_values = []
        
        for sample in samples:
            # Filter for current family, sample, and presence=1
            sample_family_data = family_data[(family_data['sample'] == sample) & (family_data['presence'] == 1)]
            present_count = len(sample_family_data)
            
            if use_proportions:
                # Calculate proportion relative to total present introners in this sample
                total_present = sample_totals[sample]
                value = present_count / total_present if total_present > 0 else 0
            else:
                # Use absolute count
                value = present_count
                
            sample_values.append(value)
            
        result_data[family] = sample_values
    
    # Create a DataFrame for easier plotting
    plot_data = pd.DataFrame(result_data, index=samples)
    
    # Set up the plot with appropriate size based on number of families
    plt.figure(figsize=(15, 10))
    
    # Set up positions for grouped bars
    x = np.arange(len(top_families))
    width = 0.8 / len(samples)  # Width of each bar, adjusted for number of samples
    
    # Plot bars for each sample
    for i, sample in enumerate(samples):
        positions = x + (i - len(samples)/2 + 0.5) * width
        plt.bar(positions, plot_data.loc[sample], width=width, label=sample, alpha=0.8)
    
    # Add labels and title
    plt.xlabel('Family')
    
    if use_proportions:
        plt.ylabel('Proportion of Present Introners')
        plt.title('Proportion of Present Introners by Family Across Samples')
        plt.ylim(0, 1.0)
    else:
        plt.ylabel('Number of Present Introners')
        plt.title('Count of Present Introners by Family Across Samples')
        plt.ylim(bottom=0)
    
    plt.xticks(x, top_families, rotation=45, ha='right')
    
    # Add a grid for better readability
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    # Add a legend with smaller font size and multiple columns if needed
    legend_cols = 1 if len(samples) <= 10 else 2
    plt.legend(title="Sample", fontsize=8, title_fontsize=10, ncol=legend_cols)
    
    plt.tight_layout()
    
    # Save the plot with appropriate suffix
    suffix = "_proportion" if use_proportions else "_count"
    plt.savefig(f"figures/{output_prefix}_family_distribution{suffix}.png")
    print(f"Family distribution plot saved to figures/{output_prefix}_family_distribution{suffix}.png")

def plot_statistics(df, pivot_df, group1, group2, afs_2d, output_prefix="genotype_stats"):
    """Generate plots for the genotype matrix statistics."""
    # Create a directory for plots if it doesn't exist
    import os
    if not os.path.exists("figures"):
        os.makedirs("figures")
    
    # Presence distribution
    plt.figure(figsize=(10, 6))
    presence_counts = df['presence'].value_counts().sort_index()
    sns.barplot(x=presence_counts.index, y=presence_counts.values)
    plt.xlabel('Presence Value')
    plt.ylabel('Count')
    plt.title('Distribution of Presence Values')
    plt.xticks([0, 1, 2], ['1 (Present)', '2 (Absent)', '3 (Missing)'])
    plt.tight_layout()
    plt.savefig(f"figures/{output_prefix}_presence_distribution.png")
    
    # Create a population column
    df['population'] = df['sample'].apply(lambda x: 'Group 2' if x in group2 else 'Group 1')
    
    # Presence distribution by population
    plt.figure(figsize=(10, 6))
    pop_presence = pd.crosstab(df['population'], df['presence'])
    pop_presence.plot(kind='bar', stacked=False)
    plt.xlabel('Population')
    plt.ylabel('Count')
    plt.title('Presence Values by Population')
    plt.legend(['Present', 'Absent', 'Missing'])
    plt.tight_layout()
    plt.savefig(f"figures/{output_prefix}_presence_by_population.png")
    
    # Heatmap of presence/absence
    if len(pivot_df) <= 100:  # Only generate heatmap for small matrices
        plt.figure(figsize=(12, 10))
        sns.heatmap(pivot_df, cmap='viridis', cbar_kws={'label': 'Presence Value'})
        plt.title('Genotype Matrix Heatmap')
        plt.tight_layout()
        plt.savefig(f"figures/{output_prefix}_heatmap.png")
    
    # Plot the family distributions - both count and proportion versions
    plot_introner_family_distribution(df, output_prefix, use_proportions=False)  # Absolute counts
    plot_introner_family_distribution(df, output_prefix, use_proportions=True)   # Proportions
    
    # 2D Allele Frequency Spectrum
    plt.figure(figsize=(10, 8))
    
    # Plot as a heatmap with log scale for better visualization
    # Add 1 to avoid log(0)
    with np.errstate(divide='ignore'):
        log_afs = np.log10(afs_2d + 1)
    
    # Mask cells with zero
    mask = afs_2d == 0
    
    sns.heatmap(log_afs, annot=afs_2d, fmt='d', cmap='viridis', mask=mask,
                cbar_kws={'label': 'log10(count + 1)'})
    plt.xlabel('Number of Group 2 Samples with Introner Present')
    plt.ylabel('Number of Group 1 Samples with Introner Present')
    plt.title('2D Allele Frequency Spectrum')
    plt.tight_layout()
    plt.savefig(f"figures/{output_prefix}_2d_afs.png")
    
    # 1D Allele Frequency Spectrum for each group
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Group 1 AFS (sum across Group 2)
    g1_afs = np.sum(afs_2d, axis=1)
    ax1.bar(range(len(g1_afs)), g1_afs)
    ax1.set_xlabel('Number of Group 1 Samples with Introner Present')
    ax1.set_ylabel('Count')
    ax1.set_title('Group 1 Allele Frequency Spectrum')
    
    # Group 2 AFS (sum across Group 1)
    g2_afs = np.sum(afs_2d, axis=0)
    ax2.bar(range(len(g2_afs)), g2_afs)
    ax2.set_xlabel('Number of Group 2 Samples with Introner Present')
    ax2.set_ylabel('Count')
    ax2.set_title('Group 2 Allele Frequency Spectrum')
    
    plt.tight_layout()
    plt.savefig(f"figures/{output_prefix}_1d_afs.png")
    
    # Per-sample statistics
    sample_stats = pd.DataFrame({
        'Present (1)': (pivot_df == 1).sum(),
        'Absent (2)': (pivot_df == 2).sum(),
        'Missing (3)': (pivot_df == 3).sum()
    })
    sample_stats['Population'] = ['Group 2' if s in group2 else 'Group 1' for s in sample_stats.index]
    
    # Reset index for easier plotting - make sure to use 'sample' (lowercase) to match later code
    sample_stats = sample_stats.reset_index().rename(columns={'index': 'sample'})
    
    # Sort by population and then by sample name - use consistent column names
    sample_stats = sample_stats.sort_values(['Population', 'sample'])
    
    # Melt the data for easier plotting
    melted = pd.melt(sample_stats, id_vars=['sample', 'Population'], 
                    value_vars=['Present (1)', 'Absent (2)', 'Missing (3)'],
                    var_name='Status', value_name='Count')
    
    # Create the plot
    plt.figure(figsize=(14, 8))
    sns.barplot(x='sample', y='Count', hue='Status', data=melted)
    plt.xticks(rotation=90)
    plt.title('Presence/Absence/Missing by Sample and Population')
    plt.tight_layout()
    plt.savefig(f"figures/{output_prefix}_sample_population_distribution.png")
    
    print(f"\nPlots saved to figures/{output_prefix}_*.png")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python analyze_genotype_matrix.py <genotype_matrix.tsv> [output_prefix]")
        sys.exit(1)
    
    file_path = sys.argv[1]
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else "genotype_stats"
    
    df, pivot_df, group1, group2, afs_2d = analyze_genotype_matrix(file_path)
    plot_statistics(df, pivot_df, group1, group2, afs_2d, output_prefix)