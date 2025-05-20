#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser(description='Visualize diversity metrics')
    parser.add_argument('--input', help='Input file with diversity metrics')
    parser.add_argument('--output', help='Output file for visualization')
    return parser.parse_args()

def plot_monomorphic_diversity(df, output_file):
    """Create plots for monomorphic groups"""
    # Filter for monomorphic groups
    mono_df = df[df['category'] == 'monomorphic_both_groups'].copy()
    
    if mono_df.empty:
        print("No monomorphic groups to plot")
        return
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    
    # Plot 1: dXY by group status
    sns.boxplot(x='g1_status', y='dxy', hue='g2_status', data=mono_df, ax=axes[0, 0])
    axes[0, 0].set_title('dXY by Group Status')
    axes[0, 0].set_xlabel('Group 1 Status')
    axes[0, 0].set_ylabel('dXY')
    
    # Plot 2: Pi by group and status
    g1_pi = mono_df[['g1_status', 'g1_pi']].copy()
    g1_pi['group'] = 'Group 1'
    g1_pi.rename(columns={'g1_status': 'status', 'g1_pi': 'pi'}, inplace=True)
    
    g2_pi = mono_df[['g2_status', 'g2_pi']].copy()
    g2_pi['group'] = 'Group 2'
    g2_pi.rename(columns={'g2_status': 'status', 'g2_pi': 'pi'}, inplace=True)
    
    pi_df = pd.concat([g1_pi, g2_pi])
    
    sns.boxplot(x='status', y='pi', hue='group', data=pi_df, ax=axes[0, 1])
    axes[0, 1].set_title('Pi by Group and Status')
    axes[0, 1].set_xlabel('Status')
    axes[0, 1].set_ylabel('Pi')
    
    # Plot 3: Status combination counts
    status_counts = mono_df.groupby(['g1_status', 'g2_status']).size().reset_index()
    status_counts.columns = ['g1_status', 'g2_status', 'count']
    status_counts['status_combo'] = status_counts['g1_status'] + '/' + status_counts['g2_status']
    
    sns.barplot(x='status_combo', y='count', data=status_counts, ax=axes[0, 2])
    axes[0, 2].set_title('Counts by Status Combination')
    axes[0, 2].set_xlabel('Status (G1/G2)')
    axes[0, 2].set_ylabel('Count')
    axes[0, 2].set_xticklabels(axes[0, 2].get_xticklabels(), rotation=45)
    
    # Plot 4: dXY vs Pi (Group 1)
    sns.scatterplot(x='g1_pi', y='dxy', hue='g1_status', style='g2_status', data=mono_df, ax=axes[1, 0])
    axes[1, 0].set_title('dXY vs Pi (Group 1)')
    axes[1, 0].set_xlabel('Pi (Group 1)')
    axes[1, 0].set_ylabel('dXY')
    
    # Plot 5: dXY vs Pi (Group 2)
    sns.scatterplot(x='g2_pi', y='dxy', hue='g2_status', style='g1_status', data=mono_df, ax=axes[1, 1])
    axes[1, 1].set_title('dXY vs Pi (Group 2)')
    axes[1, 1].set_xlabel('Pi (Group 2)')
    axes[1, 1].set_ylabel('dXY')
    
    # Plot 6: dXY distribution by status combination
    # Create status combination column
    mono_df['status_combo'] = mono_df['g1_status'] + '/' + mono_df['g2_status']
    
    sns.boxplot(x='status_combo', y='dxy', data=mono_df, ax=axes[1, 2])
    axes[1, 2].set_title('dXY by Status Combination')
    axes[1, 2].set_xlabel('Status (G1/G2)')
    axes[1, 2].set_ylabel('dXY')
    axes[1, 2].set_xticklabels(axes[1, 2].get_xticklabels(), rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def plot_polymorphic_diversity(df, output_file):
    """Create plots for polymorphic group1"""
    # Filter for polymorphic group1
    poly_df = df[df['category'] == 'polymorphic_group1'].copy()
    
    if poly_df.empty:
        print("No polymorphic group1 to plot")
        return
    
    # Create figure with multiple subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Plot 1: Pi within present vs absent
    present_pi = poly_df[['ortholog_id', 'g1_present_pi']].copy()
    present_pi['category'] = 'Present'
    present_pi.rename(columns={'g1_present_pi': 'pi'}, inplace=True)
    
    absent_pi = poly_df[['ortholog_id', 'g1_absent_pi']].copy()
    absent_pi['category'] = 'Absent'
    absent_pi.rename(columns={'g1_absent_pi': 'pi'}, inplace=True)
    
    pi_df = pd.concat([present_pi, absent_pi])
    
    sns.boxplot(x='category', y='pi', data=pi_df, ax=axes[0, 0])
    axes[0, 0].set_title('Pi within Present vs Absent (Group 1)')
    axes[0, 0].set_xlabel('Category')
    axes[0, 0].set_ylabel('Pi')
    
    # Plot 2: Pi between present and absent
    sns.histplot(poly_df['g1_present_absent_pi'].dropna(), ax=axes[0, 1], kde=True)
    axes[0, 1].set_title('Pi between Present and Absent (Group 1)')
    axes[0, 1].set_xlabel('Pi')
    axes[0, 1].set_ylabel('Count')
    
    # Plot 3: Pi by sample count
    poly_df['total_samples'] = poly_df['g1_present_count'] + poly_df['g1_absent_count']
    sns.scatterplot(x='total_samples', y='g1_present_absent_pi', hue='g2_status', data=poly_df, ax=axes[1, 0])
    axes[1, 0].set_title('Pi between Present/Absent by Sample Count')
    axes[1, 0].set_xlabel('Total Samples in Group 1')
    axes[1, 0].set_ylabel('Pi between Present/Absent')
    
    # Plot 4: Pi by present/absent ratio
    poly_df['present_ratio'] = poly_df['g1_present_count'] / poly_df['total_samples']
    sns.scatterplot(x='present_ratio', y='g1_present_absent_pi', hue='g2_status', data=poly_df, ax=axes[1, 1])
    axes[1, 1].set_title('Pi between Present/Absent by Present Ratio')
    axes[1, 1].set_xlabel('Present Ratio in Group 1')
    axes[1, 1].set_ylabel('Pi between Present/Absent')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()

def create_status_summary_table(df):
    """Create summary table for ortholog status combinations"""
    # Create status combination
    df['status_combo'] = df['g1_status'] + '/' + df['g2_status']
    
    # Group by status combination and calculate metrics
    status_summary = df.groupby(['category', 'status_combo']).agg({
        'ortholog_id': 'count',
        'dxy': ['mean', 'std', 'min', 'max'],
        'g1_pi': ['mean', 'std'],
        'g2_pi': ['mean', 'std']
    }).reset_index()
    
    # Flatten column names
    status_summary.columns = [
        '_'.join(col).strip('_') for col in status_summary.columns.values
    ]
    
    # Rename columns for clarity
    status_summary.rename(columns={
        'ortholog_id_count': 'Count',
        'dxy_mean': 'Mean dXY',
        'dxy_std': 'Std dXY',
        'dxy_min': 'Min dXY',
        'dxy_max': 'Max dXY',
        'g1_pi_mean': 'Mean Pi (G1)',
        'g1_pi_std': 'Std Pi (G1)',
        'g2_pi_mean': 'Mean Pi (G2)',
        'g2_pi_std': 'Std Pi (G2)'
    }, inplace=True)
    
    return status_summary

def main():
    args = parse_arguments()
    
    # Load diversity metrics
    print(f"Loading diversity metrics from {args.input}")
    df = pd.read_csv(args.input, sep="\t")
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(args.output)
    os.makedirs(output_dir, exist_ok=True)
    
    # Create base filename without extension
    base_output = os.path.splitext(args.output)[0]
    
    # Plot monomorphic diversity
    mono_output = f"{base_output}_monomorphic.png"
    print(f"Creating monomorphic diversity plots: {mono_output}")
    plot_monomorphic_diversity(df, mono_output)
    
    # Plot polymorphic diversity
    poly_output = f"{base_output}_polymorphic.png"
    print(f"Creating polymorphic diversity plots: {poly_output}")
    plot_polymorphic_diversity(df, poly_output)
    
    # Create combined summary plot
    print(f"Creating combined summary plot: {args.output}")
    
    # Create figure with summary statistics
    plt.figure(figsize=(14, 10))
    
    # Count ortholog groups by category
    category_counts = df['category'].value_counts()
    
    # Create status summary table
    status_summary = create_status_summary_table(df)
    
    # Plot summary statistics
    plt.subplot(2, 2, 1)
    plt.pie(category_counts, labels=category_counts.index, autopct='%1.1f%%')
    plt.title('Ortholog Groups by Category')
    
    # Create status combination counts for monomorphic groups
    mono_df = df[df['category'] == 'monomorphic_both_groups'].copy()
    if not mono_df.empty:
        mono_df['status_combo'] = mono_df['g1_status'] + '/' + mono_df['g2_status']
        status_combo_counts = mono_df['status_combo'].value_counts()
        
        plt.subplot(2, 2, 2)
        plt.pie(status_combo_counts, labels=status_combo_counts.index, autopct='%1.1f%%')
        plt.title('Monomorphic Groups by Status Combination')
    
    # Summary statistics for dXY
    plt.subplot(2, 2, 3)
    mono_status_combos = mono_df['status_combo'].unique() if not mono_df.empty else []
    
    dxy_data = []
    dxy_labels = []
    
    # Collect dXY data for different status combinations
    for combo in mono_status_combos:
        combo_data = mono_df[mono_df['status_combo'] == combo]['dxy'].dropna()
        if not combo_data.empty:
            dxy_data.append(combo_data)
            dxy_labels.append(combo)
    
    if dxy_data:
        plt.boxplot(dxy_data, labels=dxy_labels)
        plt.title('dXY by Status Combination')
        plt.ylabel('dXY')
        plt.xticks(rotation=45)
    
    # Summary statistics table
    plt.subplot(2, 2, 4)
    plt.axis('off')
    
    # Format the summary table for display
    display_cols = ['category', 'status_combo', 'Count', 'Mean dXY', 'Mean Pi (G1)', 'Mean Pi (G2)']
    display_summary = status_summary[display_cols].copy()
    
    # Format numbers to reduce decimal places
    for col in ['Mean dXY', 'Mean Pi (G1)', 'Mean Pi (G2)']:
        display_summary[col] = display_summary[col].map(lambda x: f"{x:.4f}" if pd.notnull(x) else "NA")
    
    # Create a table
    table = plt.table(
        cellText=display_summary.values,
        colLabels=display_summary.columns,
        loc='center',
        cellLoc='center'
    )
    
    # Adjust table appearance
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.5)
    
    plt.title('Summary Statistics by Status Combination')
    
    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()
    
    # Save detailed status summary to TSV
    summary_file = f"{base_output}_status_summary.tsv"
    status_summary.to_csv(summary_file, sep="\t", index=False)
    print(f"Saved detailed status summary to {summary_file}")
    
    print("Done!")

if __name__ == "__main__":
    main()