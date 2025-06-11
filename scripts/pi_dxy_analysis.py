#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from scipy import stats
from matplotlib.patches import Rectangle

def load_and_prepare_data(filepath):
    """Load diversity metrics and prepare for analysis"""
    df = pd.read_csv(filepath, sep='\t')
    
    # Create a combined status column for easier grouping
    df['status_combination'] = df['g1_status'] + '_' + df['g2_status']
    
    # Calculate allele frequencies for polymorphic groups
    df['g1_present_freq'] = df['g1_present_count'] / (df['g1_present_count'] + df['g1_absent_count'])
    df['g1_total_samples'] = df['g1_present_count'] + df['g1_absent_count']
    
    return df

def plot_pi_dxy_overview(df, output_prefix):
    """Create comprehensive Pi vs dXY plots"""
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Introner Diversity Patterns: Pi vs dXY Analysis', fontsize=16, fontweight='bold')
    
    # Filter for monomorphic groups with valid dXY values
    mono_df = df[(df['category'] == 'monomorphic_both_groups') & 
                 (df['dxy'].notna()) & (df['dxy'] > 0)].copy()
    
    # Plot 1: G1 Pi vs dXY by status combination
    ax1 = axes[0, 0]
    if not mono_df.empty:
        for status_combo in mono_df['status_combination'].unique():
            subset = mono_df[mono_df['status_combination'] == status_combo]
            if len(subset) > 0:
                ax1.scatter(subset['g1_pi'], subset['dxy'], 
                           label=status_combo.replace('_', ' + '), 
                           alpha=0.7, s=50)
        
        ax1.set_xlabel('Pi (Group 1)')
        ax1.set_ylabel('dXY (Between Groups)')
        ax1.set_title('Group 1 Pi vs dXY by Status')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Add diagonal line for reference (Pi = dXY)
        max_val = max(mono_df['g1_pi'].max(), mono_df['dxy'].max())
        ax1.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Pi = dXY')
    
    # Plot 2: G2 Pi vs dXY by status combination
    ax2 = axes[0, 1]
    if not mono_df.empty:
        for status_combo in mono_df['status_combination'].unique():
            subset = mono_df[mono_df['status_combination'] == status_combo]
            if len(subset) > 0:
                ax2.scatter(subset['g2_pi'], subset['dxy'], 
                           label=status_combo.replace('_', ' + '), 
                           alpha=0.7, s=50)
        
        ax2.set_xlabel('Pi (Group 2)')
        ax2.set_ylabel('dXY (Between Groups)')
        ax2.set_title('Group 2 Pi vs dXY by Status')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # Add diagonal line for reference
        max_val = max(mono_df['g2_pi'].max(), mono_df['dxy'].max())
        ax2.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Pi = dXY')
    
    # Plot 3: G1 vs G2 Pi colored by dXY
    ax3 = axes[0, 2]
    if not mono_df.empty:
        scatter = ax3.scatter(mono_df['g1_pi'], mono_df['g2_pi'], 
                             c=mono_df['dxy'], cmap='viridis', 
                             alpha=0.7, s=50)
        ax3.set_xlabel('Pi (Group 1)')
        ax3.set_ylabel('Pi (Group 2)')
        ax3.set_title('Group 1 vs Group 2 Pi\n(colored by dXY)')
        
        # Add diagonal line
        max_val = max(mono_df['g1_pi'].max(), mono_df['g2_pi'].max())
        ax3.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='G1 Pi = G2 Pi')
        ax3.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax3)
        cbar.set_label('dXY', rotation=270, labelpad=15)
    
    # Plot 4: Polymorphic Group 1 analysis
    ax4 = axes[1, 0]
    poly_df = df[(df['category'] == 'polymorphic_group1') & 
                 (df['g1_present_pi'].notna()) & (df['g1_absent_pi'].notna())].copy()
    
    if not poly_df.empty:
        ax4.scatter(poly_df['g1_present_pi'], poly_df['g1_absent_pi'], 
                   c=poly_df['g1_present_absent_pi'], cmap='plasma', 
                   alpha=0.7, s=50)
        ax4.set_xlabel('Pi (Present Samples)')
        ax4.set_ylabel('Pi (Absent Samples)')
        ax4.set_title('Polymorphic Group 1:\nPresent vs Absent Pi')
        ax4.grid(True, alpha=0.3)
        
        # Add diagonal line
        max_val = max(poly_df['g1_present_pi'].max(), poly_df['g1_absent_pi'].max())
        ax4.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
        
        # Add colorbar
        scatter4 = ax4.collections[0]
        cbar4 = plt.colorbar(scatter4, ax=ax4)
        cbar4.set_label('Pi (Present vs Absent)', rotation=270, labelpad=15)
    
    # Plot 5: dXY distribution by status combination
    ax5 = axes[1, 1]
    if not mono_df.empty:
        status_order = sorted(mono_df['status_combination'].unique())
        box_data = [mono_df[mono_df['status_combination'] == status]['dxy'].dropna() 
                   for status in status_order]
        
        bp = ax5.boxplot(box_data, labels=[s.replace('_', '\n') for s in status_order], 
                        patch_artist=True)
        
        # Color the boxes
        colors = sns.color_palette("husl", len(status_order))
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax5.set_ylabel('dXY (Between Groups)')
        ax5.set_title('dXY Distribution by\nStatus Combination')
        ax5.grid(True, alpha=0.3)
    
    # Plot 6: Allele frequency vs diversity in polymorphic groups
    ax6 = axes[1, 2]
    if not poly_df.empty:
        # Plot present vs absent Pi by allele frequency
        ax6.scatter(poly_df['g1_present_freq'], poly_df['g1_present_pi'], 
                   label='Present Samples', alpha=0.7, s=50, color='blue')
        ax6.scatter(poly_df['g1_present_freq'], poly_df['g1_absent_pi'], 
                   label='Absent Samples', alpha=0.7, s=50, color='red')
        ax6.scatter(poly_df['g1_present_freq'], poly_df['g1_present_absent_pi'], 
                   label='Present vs Absent', alpha=0.7, s=50, color='green')
        
        ax6.set_xlabel('Allele Frequency (Present)')
        ax6.set_ylabel('Pi')
        ax6.set_title('Diversity vs Allele Frequency\n(Polymorphic Group 1)')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_pi_dxy_overview.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_pi_dxy_overview.pdf', bbox_inches='tight')
    plt.show()

def plot_evolutionary_signatures(df, output_prefix):
    """Plot patterns that might indicate different evolutionary processes"""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Evolutionary Signatures in Introner Diversity', fontsize=16, fontweight='bold')
    
    # Filter monomorphic groups
    mono_df = df[(df['category'] == 'monomorphic_both_groups') & 
                 (df['dxy'].notna()) & (df['dxy'] > 0)].copy()
    
    # Plot 1: Pi/dXY ratio analysis
    ax1 = axes[0, 0]
    if not mono_df.empty:
        # Calculate Pi/dXY ratios
        mono_df['g1_pi_dxy_ratio'] = mono_df['g1_pi'] / mono_df['dxy']
        mono_df['g2_pi_dxy_ratio'] = mono_df['g2_pi'] / mono_df['dxy']
        
        for status_combo in mono_df['status_combination'].unique():
            subset = mono_df[mono_df['status_combination'] == status_combo]
            if len(subset) > 0:
                ax1.scatter(subset['g1_pi_dxy_ratio'], subset['g2_pi_dxy_ratio'], 
                           label=status_combo.replace('_', ' + '), alpha=0.7, s=50)
        
        ax1.set_xlabel('Pi/dXY Ratio (Group 1)')
        ax1.set_ylabel('Pi/dXY Ratio (Group 2)')
        ax1.set_title('Pi/dXY Ratios\n(Low ratios suggest recent selection)')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Add reference lines
        ax1.axhline(y=1, color='k', linestyle='--', alpha=0.5)
        ax1.axvline(x=1, color='k', linestyle='--', alpha=0.5)
    
    # Plot 2: Neutrality test approximation
    ax2 = axes[0, 1]
    if not mono_df.empty:
        # Rough neutrality expectation: Pi should be proportional to dXY
        for status_combo in mono_df['status_combination'].unique():
            subset = mono_df[mono_df['status_combination'] == status_combo]
            if len(subset) > 0:
                # Calculate deviation from neutral expectation
                expected_pi = subset['dxy'] * 0.5  # Rough approximation
                deviation_g1 = (subset['g1_pi'] - expected_pi) / expected_pi
                deviation_g2 = (subset['g2_pi'] - expected_pi) / expected_pi
                
                ax2.scatter(deviation_g1, deviation_g2, 
                           label=status_combo.replace('_', ' + '), alpha=0.7, s=50)
        
        ax2.set_xlabel('Deviation from Neutral (Group 1)')
        ax2.set_ylabel('Deviation from Neutral (Group 2)')
        ax2.set_title('Deviation from Neutral Expectation\n(Negative = reduced diversity)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        ax2.axhline(y=0, color='k', linestyle='-', alpha=0.5)
        ax2.axvline(x=0, color='k', linestyle='-', alpha=0.5)
    
    # Plot 3: Population differentiation (Fst approximation)
    ax3 = axes[1, 0]
    if not mono_df.empty:
        # Calculate approximate Fst: (dXY - average_Pi) / dXY
        mono_df['avg_pi'] = (mono_df['g1_pi'] + mono_df['g2_pi']) / 2
        mono_df['fst_approx'] = (mono_df['dxy'] - mono_df['avg_pi']) / mono_df['dxy']
        
        # Plot Fst vs average Pi
        scatter = ax3.scatter(mono_df['avg_pi'], mono_df['fst_approx'], 
                             c=[hash(x) for x in mono_df['status_combination']], 
                             cmap='tab10', alpha=0.7, s=50)
        
        ax3.set_xlabel('Average Pi (Both Groups)')
        ax3.set_ylabel('Approximate Fst')
        ax3.set_title('Population Differentiation\nvs Average Diversity')
        ax3.grid(True, alpha=0.3)
        ax3.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Plot 4: Polymorphic group evolutionary patterns
    ax4 = axes[1, 1]
    poly_df = df[(df['category'] == 'polymorphic_group1') & 
                 (df['g1_present_absent_pi'].notna())].copy()
    
    if not poly_df.empty:
        # Compare within-group diversity to between-state diversity
        poly_df['avg_within_pi'] = (poly_df['g1_present_pi'] + poly_df['g1_absent_pi']) / 2
        
        ax4.scatter(poly_df['avg_within_pi'], poly_df['g1_present_absent_pi'], 
                   alpha=0.7, s=50, color='purple')
        
        ax4.set_xlabel('Average Pi Within States')
        ax4.set_ylabel('Pi Between Present/Absent')
        ax4.set_title('Within vs Between State Diversity\n(Polymorphic Group 1)')
        ax4.grid(True, alpha=0.3)
        
        # Add diagonal line
        max_val = max(poly_df['avg_within_pi'].max(), poly_df['g1_present_absent_pi'].max())
        ax4.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_evolutionary_signatures.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_evolutionary_signatures.pdf', bbox_inches='tight')
    plt.show()

def generate_summary_statistics(df, output_prefix):
    """Generate summary statistics and save to file"""
    
    summary_stats = []
    
    # Monomorphic groups statistics
    mono_df = df[(df['category'] == 'monomorphic_both_groups') & 
                 (df['dxy'].notna()) & (df['dxy'] > 0)]
    
    if not mono_df.empty:
        for status_combo in mono_df['status_combination'].unique():
            subset = mono_df[mono_df['status_combination'] == status_combo]
            
            stats_dict = {
                'category': 'monomorphic_both_groups',
                'status_combination': status_combo,
                'n_orthologs': len(subset),
                'mean_g1_pi': subset['g1_pi'].mean(),
                'mean_g2_pi': subset['g2_pi'].mean(),
                'mean_dxy': subset['dxy'].mean(),
                'median_g1_pi': subset['g1_pi'].median(),
                'median_g2_pi': subset['g2_pi'].median(),
                'median_dxy': subset['dxy'].median(),
                'std_g1_pi': subset['g1_pi'].std(),
                'std_g2_pi': subset['g2_pi'].std(),
                'std_dxy': subset['dxy'].std()
            }
            
            # Add Pi/dXY ratios
            if subset['dxy'].mean() > 0:
                stats_dict['mean_g1_pi_dxy_ratio'] = subset['g1_pi'].mean() / subset['dxy'].mean()
                stats_dict['mean_g2_pi_dxy_ratio'] = subset['g2_pi'].mean() / subset['dxy'].mean()
            
            summary_stats.append(stats_dict)
    
    # Polymorphic groups statistics
    poly_df = df[(df['category'] == 'polymorphic_group1') & 
                 (df['g1_present_absent_pi'].notna())]
    
    if not poly_df.empty:
        for g2_status in poly_df['g2_status'].unique():
            subset = poly_df[poly_df['g2_status'] == g2_status]
            
            stats_dict = {
                'category': 'polymorphic_group1',
                'status_combination': f'polymorphic_{g2_status}',
                'n_orthologs': len(subset),
                'mean_g1_present_pi': subset['g1_present_pi'].mean(),
                'mean_g1_absent_pi': subset['g1_absent_pi'].mean(),
                'mean_g1_present_absent_pi': subset['g1_present_absent_pi'].mean(),
                'median_g1_present_pi': subset['g1_present_pi'].median(),
                'median_g1_absent_pi': subset['g1_absent_pi'].median(),
                'median_g1_present_absent_pi': subset['g1_present_absent_pi'].median(),
                'mean_allele_freq': subset['g1_present_freq'].mean(),
                'median_allele_freq': subset['g1_present_freq'].median()
            }
            
            summary_stats.append(stats_dict)
    
    # Save to file
    summary_df = pd.DataFrame(summary_stats)
    summary_df.to_csv(f'{output_prefix}_summary_statistics.tsv', sep='\t', index=False)
    
    # Print key findings
    print("\n=== DIVERSITY ANALYSIS SUMMARY ===")
    print(f"Total ortholog groups analyzed: {len(df)}")
    print(f"Monomorphic groups: {len(mono_df)}")
    print(f"Polymorphic Group 1: {len(poly_df)}")
    
    if not mono_df.empty:
        print(f"\nMonomorphic groups - Average dXY: {mono_df['dxy'].mean():.4f}")
        print(f"Monomorphic groups - Average Pi (G1): {mono_df['g1_pi'].mean():.4f}")
        print(f"Monomorphic groups - Average Pi (G2): {mono_df['g2_pi'].mean():.4f}")
    
    if not poly_df.empty:
        print(f"\nPolymorphic Group 1 - Average Pi (present): {poly_df['g1_present_pi'].mean():.4f}")
        print(f"Polymorphic Group 1 - Average Pi (absent): {poly_df['g1_absent_pi'].mean():.4f}")
        print(f"Polymorphic Group 1 - Average Pi (present vs absent): {poly_df['g1_present_absent_pi'].mean():.4f}")
    
    return summary_df

def main():
    parser = argparse.ArgumentParser(description='Analyze Pi vs dXY patterns in introner diversity data')
    parser.add_argument('input_file', help='Input TSV file with diversity metrics')
    parser.add_argument('--output_prefix', default='introner_diversity', 
                       help='Prefix for output files (default: introner_diversity)')
    
    args = parser.parse_args()
    
    # Load and prepare data
    print("Loading diversity data...")
    df = load_and_prepare_data(args.input_file)
    
    # Generate plots
    print("Creating Pi vs dXY overview plots...")
    plot_pi_dxy_overview(df, args.output_prefix)
    
    print("Creating evolutionary signature plots...")
    plot_evolutionary_signatures(df, args.output_prefix)
    
    # Generate summary statistics
    print("Generating summary statistics...")
    summary_df = generate_summary_statistics(df, args.output_prefix)
    
    print(f"\nAnalysis complete! Output files saved with prefix: {args.output_prefix}")
    print("Generated files:")
    print(f"  - {args.output_prefix}_pi_dxy_overview.png/pdf")
    print(f"  - {args.output_prefix}_evolutionary_signatures.png/pdf") 
    print(f"  - {args.output_prefix}_summary_statistics.tsv")

if __name__ == "__main__":
    main()