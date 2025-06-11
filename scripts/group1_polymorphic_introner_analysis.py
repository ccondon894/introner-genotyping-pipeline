#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
from scipy import stats
from scipy.stats import pearsonr, spearmanr
import warnings
warnings.filterwarnings('ignore')

def load_and_filter_data(filepath):
    """Load diversity data and filter for polymorphic Group 1"""
    df = pd.read_csv(filepath, sep='\t')
    
    # Filter for polymorphic Group 1 with valid data
    poly_df = df[
        (df['category'] == 'polymorphic_group1') & 
        (df['g1_present_pi'].notna()) & 
        (df['g1_absent_pi'].notna()) & 
        (df['g1_present_absent_pi'].notna()) &
        (df['g1_present_count'] >= 2) &  # At least 2 present
        (df['g1_absent_count'] >= 2) &   # At least 2 absent
        (df['g1_present_count'] <= 9)    # Max 9 present (out of 11 total)
    ].copy()
    
    # Calculate derived metrics
    poly_df['total_samples'] = poly_df['g1_present_count'] + poly_df['g1_absent_count']
    poly_df['present_freq'] = poly_df['g1_present_count'] / poly_df['total_samples']
    poly_df['pi_ratio_between_within'] = poly_df['g1_present_absent_pi'] / ((poly_df['g1_present_pi'] + poly_df['g1_absent_pi']) / 2)
    
    return poly_df

def calculate_afs_statistics(df):
    """Calculate statistics for each allele frequency class"""
    afs_stats = []
    
    for present_count in range(2, 10):  # 2 to 9 present alleles
        subset = df[df['g1_present_count'] == present_count]
        
        if len(subset) == 0:
            continue
            
        stats_dict = {
            'present_count': present_count,
            'absent_count': subset['g1_absent_count'].iloc[0] if len(subset) > 0 else None,
            'present_freq': present_count / 11,  # Assuming 11 total samples in G1
            'n_loci': len(subset),
            
            # Basic Pi statistics
            'mean_pi_present': subset['g1_present_pi'].mean(),
            'mean_pi_absent': subset['g1_absent_pi'].mean(),
            'mean_pi_between': subset['g1_present_absent_pi'].mean(),
            
            'median_pi_present': subset['g1_present_pi'].median(),
            'median_pi_absent': subset['g1_absent_pi'].median(),
            'median_pi_between': subset['g1_present_absent_pi'].median(),
            
            'std_pi_present': subset['g1_present_pi'].std(),
            'std_pi_absent': subset['g1_absent_pi'].std(),
            'std_pi_between': subset['g1_present_absent_pi'].std(),
            
            # Ratios and comparisons
            'mean_pi_ratio': subset['pi_ratio_between_within'].mean(),
            'median_pi_ratio': subset['pi_ratio_between_within'].median(),
            
            # Statistical tests
            'pi_present_vs_absent_pvalue': stats.wilcoxon(subset['g1_present_pi'], subset['g1_absent_pi'])[1] if len(subset) > 5 else np.nan,
            'pi_between_vs_within_effect_size': (subset['g1_present_absent_pi'] - (subset['g1_present_pi'] + subset['g1_absent_pi'])/2).mean(),
        }
        
        afs_stats.append(stats_dict)
    
    return pd.DataFrame(afs_stats)

def plot_afs_pi_patterns(df, afs_stats, output_prefix):
    """Create comprehensive plots of Pi patterns across allele frequencies"""
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create main figure
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle('Pi Patterns Across Present Count Spectrum\n(Polymorphic Group 1)', 
                 fontsize=16, fontweight='bold')
    
    # Plot 1: Pi by present count - line plot
    ax1 = axes[0, 0]
    if not afs_stats.empty:
        ax1.plot(afs_stats['present_count'], afs_stats['mean_pi_present'], 
                'o-', label='Pi (Present)', linewidth=2, markersize=6)
        ax1.plot(afs_stats['present_count'], afs_stats['mean_pi_absent'], 
                'o-', label='Pi (Absent)', linewidth=2, markersize=6)
        ax1.plot(afs_stats['present_count'], afs_stats['mean_pi_between'], 
                'o-', label='Pi (Present vs Absent)', linewidth=2, markersize=6)
        
        ax1.set_xlabel('G1 Present Count')
        ax1.set_ylabel('Mean Pi')
        ax1.set_title('Pi by Present Count')
        ax1.set_xticks(afs_stats['present_count'])
        ax1.legend()
        ax1.grid(True, alpha=0.3)
    
    # Plot 2: Pi ratio by present count
    ax2 = axes[0, 1]
    if not afs_stats.empty:
        ax2.plot(afs_stats['present_count'], afs_stats['mean_pi_ratio'], 
                'o-', color='purple', linewidth=2, markersize=6)
        ax2.axhline(y=1, color='red', linestyle='--', alpha=0.7, label='Equal diversity')
        
        ax2.set_xlabel('G1 Present Count')
        ax2.set_ylabel('Pi(Between) / Pi(Within) Ratio')
        ax2.set_title('Between vs Within State\nDiversity Ratio')
        ax2.set_xticks(afs_stats['present_count'])
        ax2.legend()
        ax2.grid(True, alpha=0.3)
    
    # Plot 3: Number of loci by present count
    ax3 = axes[0, 2]
    if not afs_stats.empty:
        bars = ax3.bar(afs_stats['present_count'], afs_stats['n_loci'], 
                      alpha=0.7, edgecolor='black')
        ax3.set_xlabel('G1 Present Count')
        ax3.set_ylabel('Number of Loci')
        ax3.set_title('Loci Count by\nPresent Count')
        ax3.set_xticks(afs_stats['present_count'])
        ax3.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, value in zip(bars, afs_stats['n_loci']):
            ax3.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5, 
                    str(int(value)), ha='center', va='bottom')
    
    # Plot 4: Scatter plot - Present vs Absent Pi by present count
    ax4 = axes[1, 0]
    colors = plt.cm.tab10(np.linspace(0, 1, 8))  # 8 colors for counts 2-9
    for i, present_count in enumerate(range(2, 10)):
        subset = df[df['g1_present_count'] == present_count]
        if len(subset) > 0:
            ax4.scatter(subset['g1_present_pi'], subset['g1_absent_pi'], 
                       label=f'Count={present_count}', alpha=0.6, s=30, color=colors[i])
    
    ax4.set_xlabel('Pi (Present)')
    ax4.set_ylabel('Pi (Absent)')
    ax4.set_title('Present vs Absent Pi\n(colored by present count)')
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    # Add diagonal line
    max_val = max(df['g1_present_pi'].max(), df['g1_absent_pi'].max())
    ax4.plot([0, max_val], [0, max_val], 'k--', alpha=0.5)
    
    # Plot 5: Between vs Within Pi scatter
    ax5 = axes[1, 1]
    df['avg_within_pi'] = (df['g1_present_pi'] + df['g1_absent_pi']) / 2
    
    for i, present_count in enumerate(range(2, 10)):
        subset = df[df['g1_present_count'] == present_count]
        if len(subset) > 0:
            ax5.scatter(subset['avg_within_pi'], subset['g1_present_absent_pi'], 
                       label=f'Count={present_count}', alpha=0.6, s=30, color=colors[i])
    
    ax5.set_xlabel('Average Within-State Pi')
    ax5.set_ylabel('Between-State Pi')
    ax5.set_title('Between vs Within State Pi\n(colored by present count)')
    ax5.grid(True, alpha=0.3)
    
    # Add diagonal line
    max_val = max(df['avg_within_pi'].max(), df['g1_present_absent_pi'].max())
    ax5.plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='Equal diversity')
    
    # Plot 6: Box plot of Pi ratios by present count
    ax6 = axes[1, 2]
    present_counts = sorted(df['g1_present_count'].unique())
    box_data = []
    box_labels = []
    
    for present_count in present_counts:
        subset = df[df['g1_present_count'] == present_count]
        if len(subset) > 0:
            box_data.append(subset['pi_ratio_between_within'])
            box_labels.append(str(present_count))
    
    if box_data:
        bp = ax6.boxplot(box_data, labels=box_labels, patch_artist=True)
        
        # Color the boxes with different colors
        colors_box = plt.cm.viridis(np.linspace(0, 1, len(box_data)))
        for patch, color in zip(bp['boxes'], colors_box):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        ax6.set_xlabel('G1 Present Count')
        ax6.set_ylabel('Pi(Between) / Pi(Within) Ratio')
        ax6.set_title('Diversity Ratio Distribution\nby Present Count')
        ax6.axhline(y=1, color='red', linestyle='--', alpha=0.7)
        ax6.grid(True, alpha=0.3)
    
    # Plot 7: Expected vs Observed diversity under neutrality
    ax7 = axes[2, 0]
    if not afs_stats.empty:
        # Calculate expected heterozygosity under neutrality
        afs_stats['expected_het'] = 2 * afs_stats['present_freq'] * (1 - afs_stats['present_freq'])
        
        ax7.plot(afs_stats['present_count'], afs_stats['expected_het'], 
                'o-', label='Expected (Neutral)', linewidth=2, markersize=6, color='gray')
        ax7.plot(afs_stats['present_count'], afs_stats['mean_pi_between'], 
                'o-', label='Observed (Between)', linewidth=2, markersize=6, color='blue')
        
        ax7.set_xlabel('G1 Present Count')
        ax7.set_ylabel('Diversity')
        ax7.set_title('Expected vs Observed\nDiversity')
        ax7.set_xticks(afs_stats['present_count'])
        ax7.legend()
        ax7.grid(True, alpha=0.3)
    
    # Plot 8: Correlation analysis - Present count vs between-state Pi
    ax8 = axes[2, 1]
    if len(df) > 10:
        # Calculate correlations
        corr_present_count = pearsonr(df['g1_present_count'], df['g1_present_absent_pi'])
        
        # Create scatter plot with jitter for better visualization
        jittered_counts = df['g1_present_count'] + np.random.normal(0, 0.1, len(df))
        scatter = ax8.scatter(jittered_counts, df['g1_present_absent_pi'], 
                             c=df['total_samples'], cmap='viridis', alpha=0.6, s=30)
        
        ax8.set_xlabel('G1 Present Count (jittered)')
        ax8.set_ylabel('Pi (Present vs Absent)')
        ax8.set_title(f'Present Count vs Between-State Pi\n(r={corr_present_count[0]:.3f}, p={corr_present_count[1]:.3f})')
        ax8.set_xticks(range(2, 10))
        ax8.grid(True, alpha=0.3)
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax8)
        cbar.set_label('Total Samples', rotation=270, labelpad=15)
    
    # Plot 9: Effect size by present count
    ax9 = axes[2, 2]
    if not afs_stats.empty:
        ax9.plot(afs_stats['present_count'], afs_stats['pi_between_vs_within_effect_size'], 
                'o-', linewidth=2, markersize=6, color='red')
        ax9.axhline(y=0, color='black', linestyle='-', alpha=0.5)
        
        ax9.set_xlabel('G1 Present Count')
        ax9.set_ylabel('Effect Size\n(Between - Within Pi)')
        ax9.set_title('Between vs Within\nDiversity Effect Size')
        ax9.set_xticks(afs_stats['present_count'])
        ax9.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_afs_pi_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(f'{output_prefix}_afs_pi_analysis.pdf', bbox_inches='tight')
    plt.show()

def analyze_frequency_dependent_patterns(df, afs_stats):
    """Analyze patterns that might be frequency-dependent"""
    
    print("\n=== ALLELE FREQUENCY SPECTRUM ANALYSIS ===")
    print(f"Total polymorphic loci analyzed: {len(df)}")
    print(f"Present count range: {df['g1_present_count'].min()} - {df['g1_present_count'].max()}")
    print(f"Sample size range: {df['total_samples'].min()} - {df['total_samples'].max()}")
    
    print("\n--- Pi Statistics by Present Count ---")
    for _, row in afs_stats.iterrows():
        print(f"Present count {int(row['present_count'])} (n={int(row['n_loci'])}):")
        print(f"  Pi present: {row['mean_pi_present']:.4f} ± {row['std_pi_present']:.4f}")
        print(f"  Pi absent:  {row['mean_pi_absent']:.4f} ± {row['std_pi_absent']:.4f}")
        print(f"  Pi between: {row['mean_pi_between']:.4f} ± {row['std_pi_between']:.4f}")
        print(f"  Ratio (between/within): {row['mean_pi_ratio']:.2f}")
        print()
    
    # Test for frequency-dependent patterns
    print("--- Frequency-Dependent Pattern Tests ---")
    
    # Correlation tests
    if not afs_stats.empty and len(afs_stats) > 3:
        # Test correlation between present count and Pi values
        corr_count_present, p_count_present = pearsonr(afs_stats['present_count'], afs_stats['mean_pi_present'])
        corr_count_absent, p_count_absent = pearsonr(afs_stats['present_count'], afs_stats['mean_pi_absent'])
        corr_count_between, p_count_between = pearsonr(afs_stats['present_count'], afs_stats['mean_pi_between'])
        corr_count_ratio, p_count_ratio = pearsonr(afs_stats['present_count'], afs_stats['mean_pi_ratio'])
        
        print(f"Present count vs Pi(present): r={corr_count_present:.3f}, p={p_count_present:.3f}")
        print(f"Present count vs Pi(absent):  r={corr_count_absent:.3f}, p={p_count_absent:.3f}")
        print(f"Present count vs Pi(between): r={corr_count_between:.3f}, p={p_count_between:.3f}")
        print(f"Present count vs Pi ratio:    r={corr_count_ratio:.3f}, p={p_count_ratio:.3f}")
    
    # Identify extreme frequencies
    print("\n--- Extreme Frequency Classes ---")
    if not afs_stats.empty:
        # Find frequency classes with highest/lowest diversity ratios
        max_ratio_idx = afs_stats['mean_pi_ratio'].idxmax()
        min_ratio_idx = afs_stats['mean_pi_ratio'].idxmin()
        
        print(f"Highest Pi ratio: Present count {int(afs_stats.loc[max_ratio_idx, 'present_count'])}, ratio = {afs_stats.loc[max_ratio_idx, 'mean_pi_ratio']:.2f}")
        print(f"Lowest Pi ratio:  Present count {int(afs_stats.loc[min_ratio_idx, 'present_count'])}, ratio = {afs_stats.loc[min_ratio_idx, 'mean_pi_ratio']:.2f}")
    
    # Test for U-shaped or other patterns
    if len(afs_stats) >= 5:
        print("\n--- Pattern Detection ---")
        # Test for quadratic relationship (U-shape or inverted-U)
        from numpy.polynomial import Polynomial
        
        # Fit quadratic to Pi ratio vs present count
        count_vals = afs_stats['present_count'].values
        ratio_vals = afs_stats['mean_pi_ratio'].values
        
        # Fit linear and quadratic models
        linear_fit = np.polyfit(count_vals, ratio_vals, 1)
        quad_fit = np.polyfit(count_vals, ratio_vals, 2)
        
        # Calculate R-squared for both
        linear_pred = np.polyval(linear_fit, count_vals)
        quad_pred = np.polyval(quad_fit, count_vals)
        
        ss_res_linear = np.sum((ratio_vals - linear_pred) ** 2)
        ss_res_quad = np.sum((ratio_vals - quad_pred) ** 2)
        ss_tot = np.sum((ratio_vals - np.mean(ratio_vals)) ** 2)
        
        r2_linear = 1 - (ss_res_linear / ss_tot)
        r2_quad = 1 - (ss_res_quad / ss_tot)
        
        print(f"Linear fit R²: {r2_linear:.3f}")
        print(f"Quadratic fit R²: {r2_quad:.3f}")
        
        if r2_quad > r2_linear + 0.1:
            print("Evidence for non-linear (quadratic) pattern in diversity ratio!")
        else:
            print("Pattern appears primarily linear.")

def main():
    parser = argparse.ArgumentParser(description='Analyze Pi patterns across allele frequency spectrum')
    parser.add_argument('input_file', help='Input TSV file with diversity metrics')
    parser.add_argument('--output_prefix', default='afs_pi_analysis', 
                       help='Prefix for output files (default: afs_pi_analysis)')
    
    args = parser.parse_args()
    
    # Load and filter data
    print("Loading and filtering polymorphic Group 1 data...")
    df = load_and_filter_data(args.input_file)
    
    if df.empty:
        print("No suitable data found for analysis!")
        return
    
    print(f"Found {len(df)} polymorphic loci for analysis")
    
    # Calculate AFS statistics
    print("Calculating allele frequency spectrum statistics...")
    afs_stats = calculate_afs_statistics(df)
    
    # Generate plots
    print("Creating allele frequency spectrum plots...")
    plot_afs_pi_patterns(df, afs_stats, args.output_prefix)
    
    # Analyze patterns
    analyze_frequency_dependent_patterns(df, afs_stats)
    
    # Save detailed statistics
    afs_stats.to_csv(f'{args.output_prefix}_afs_statistics.tsv', sep='\t', index=False)
    
    print(f"\nAnalysis complete! Output files saved with prefix: {args.output_prefix}")
    print("Generated files:")
    print(f"  - {args.output_prefix}_afs_pi_analysis.png/pdf")
    print(f"  - {args.output_prefix}_afs_statistics.tsv")

if __name__ == "__main__":
    main()