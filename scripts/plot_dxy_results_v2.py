import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys

tsv = sys.argv[1]
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
plt.savefig('figures/dxy_visualization_by_presence.png')

# If you'd like to use log transformation with a larger epsilon
# Uncomment and run this additional code
"""
# Create a second figure with log-transformed values (using log1p)
fig_log, axes_log = plt.subplots(3, 2, figsize=(15, 18))

for i, combo in enumerate(combinations):
    # Filter data for this combination
    subset = df[combo['filter']].copy()
    
    # Apply log1p transformation to avoid negative values
    subset['log_g1'] = np.log10(subset['g1'] + 1)
    subset['log_g2'] = np.log10(subset['g2'] + 1)
    subset['log_g1_vs_g2'] = np.log10(subset['g1_vs_g2'] + 1)
    
    # Get the axes for this row
    ax1, ax2 = axes_log[i]
    
    # Scatter plot
    ax1.scatter(subset['log_g1'], subset['log_g2'], alpha=0.5)
    ax1.set_xlabel('Log10(Pi_g1 + 1)')
    ax1.set_ylabel('Log10(Pi_g2 + 1)')
    max_val = max(subset['log_g1'].max() if not subset.empty else 0, 
                 subset['log_g2'].max() if not subset.empty else 0)
    ax1.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)  # diagonal line
    ax1.set_title(f'Log-transformed Pi g1 vs. g2 ({combo["name"]})')
    
    # Boxplot for this combination
    if not subset.empty:
        subset_melted = pd.melt(subset, 
                              id_vars=['ortholog_id'], 
                              value_vars=['log_g1', 'log_g2', 'log_g1_vs_g2'],
                              var_name='Group', value_name='Log10(Pi/dXY + 1)')
        sns.boxplot(data=subset_melted, x='Group', y='Log10(Pi/dXY + 1)', ax=ax2)
    ax2.set_title(f'Distribution of Log-transformed Pi/dXY Values ({combo["name"]})')

plt.tight_layout()
plt.savefig('figures/dxy_visualization_by_presence_log.png')
"""