# File: scripts/plot_pi_results.py

import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def main():
    tsv = sys.argv[1]
    png = sys.argv[2]
    
    # Read the data
    df = pd.read_csv(tsv, sep='\t', header=0)
    
    # Create a figure with multiple subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # Scatter plot
    ax1.scatter(df['pi_present'], df['pi_absent'], alpha=0.5)
    ax1.set_xlabel('Pi (Present)')
    ax1.set_ylabel('Pi (Absent)')
    max_val = max(df['pi_present'].max(), df['pi_absent'].max())
    ax1.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)  # diagonal line
    ax1.set_title('Pi Present vs Absent')
    
    # Boxplot for all three measurements
    df_melted = pd.melt(df, 
                        id_vars=['ortholog_id'], 
                        value_vars=['pi_present', 'pi_absent', 'pi_between'],
                        var_name='Metric', value_name='Pi')
    sns.boxplot(data=df_melted, x='Metric', y='Pi', ax=ax2)
    ax2.set_title('Distribution of Pi Values')
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(png)
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(df[['pi_present', 'pi_absent', 'pi_between']].describe())

if __name__ == "__main__":
    main()