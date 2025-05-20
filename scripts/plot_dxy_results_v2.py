# File: scripts/plot_dxy_results_v2.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

def main():
    tsv = sys.argv[1]
    png = sys.argv[2]
    
    # Read the data
    df = pd.read_csv(tsv, sep='\t', header=0)
    
    # Create a figure with 3 rows of subplots (one for each combination)
    fig, axes = plt.subplots(3, 2, figsize=(15, 18))
    
    # Define the three combinations to filter and plot
    combinations = [
        {'name': 'g1=present, g2=present', 'filter': (df['presence_g1'] == 'present') & (df['presence_g2'] == 'present')},
        {'name': 'g1=present, g2=absent', 'filter': (df['presence_g1'] == 'present') & (df['presence_g2'] == 'absent')},
        {'name': 'g1=absent, g2=present', 'filter': (df['presence_g1'] == 'absent') & (df['presence_g2'] == 'present')}
    ]
    
    # Loop through each combination and create plots in each row
    for i, combo in enumerate(combinations):
        # Filter data for this combination
        subset = df[combo['filter']]
        
        # Get the axes for this row
        ax1, ax2 = axes[i]
        
        # Scatter plot
        ax1.scatter(subset['g1'], subset['g2'], alpha=0.5)
        ax1.set_xlabel('Pi_g1')
        ax1.set_ylabel('Pi_g2')
        max_val = max(subset['g1'].max() if not subset.empty else 0, 
                     subset['g2'].max() if not subset.empty else 0)
        ax1.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)  # diagonal line
        ax1.set_title(f'Pi g1 vs. g2 ({combo["name"]})')
        
        # Boxplot for this combination
        if not subset.empty:
            subset_melted = pd.melt(subset, 
                                  id_vars=['ortholog_id'], 
                                  value_vars=['g1', 'g2', 'g1_vs_g2'],
                                  var_name='Group', value_name='Pi/dXY')
            sns.boxplot(data=subset_melted, x='Group', y='Pi/dXY', ax=ax2)
        ax2.set_title(f'Distribution of Pi/dXY Values ({combo["name"]})')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(png)

if __name__ == "__main__":
    main()