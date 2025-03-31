import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
tsv = sys.argv[1]

# Read the data
df = pd.read_csv(tsv, sep='\t', header=0)

# epsilon = 1
# df['g1'] = np.log10(df['g1'] + epsilon)
# df['g2'] = np.log10(df['g2'] + epsilon)
# df['g1_vs_g2'] = np.log10(df['g1_vs_g2'] + epsilon)

# Create a figure with multiple subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))


# Scatter plot
ax1.scatter(df['g1'], df['g2'], alpha=0.5)
ax1.set_xlabel('Pi_g1')
ax1.set_ylabel('Pi_g2')
max_val = max(df['g1'].max(), df['g2'].max())
ax1.plot([0, max_val], [0, max_val], 'r--', alpha=0.5)  # diagonal line
ax1.set_title('Pi g1 vs. g2')

# Boxplot for all three measurements
df_melted = pd.melt(df, 
                    id_vars=['ortholog_id'], 
                    value_vars=['g1', 'g2', 'g1_vs_g2'],
                    var_name='Group', value_name='Pi/dXY')
sns.boxplot(data=df_melted, x='Group', y='Pi/dXY', ax=ax2)
ax2.set_title('Distribution of Pi/dXY Values')

# Adjust layout and save
plt.tight_layout()
plt.savefig('figures/dxy_visualization.png')

# Print summary statistics
print("\nSummary Statistics:")
print(df[['g1', 'g2', 'g1_vs_g2']].describe())



# # Scatter plot with log values
# ax1.scatter(df['log_pi_present'], df['log_pi_absent'], alpha=0.5)
# ax1.set_xlabel('Log10 Pi (Present)')
# ax1.set_ylabel('Log10 Pi (Absent)')
# max_val = max(df['log_pi_present'].max(), df['log_pi_absent'].max())
# min_val = min(df['log_pi_present'].min(), df['log_pi_absent'].min())
# ax1.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)  # diagonal line
# ax1.set_title('Log10 Pi Present vs Absent')

# # Boxplot for all three log measurements
# df_melted = pd.melt(df, 
#                     id_vars=['ortholog_id'], 
#                     value_vars=['log_pi_present', 'log_pi_absent', 'log_pi_between'],
#                     var_name='Metric', value_name='Log10 Pi')

# # Clean up metric names for plotting
# df_melted['Metric'] = df_melted['Metric'].map({
#     'log_pi_present': 'Present',
#     'log_pi_absent': 'Absent',
#     'log_pi_between': 'Between'
# })

# sns.boxplot(data=df_melted, x='Metric', y='Log10 Pi', ax=ax2)
# ax2.set_title('Distribution of Log10 Pi Values')


# # Adjust layout and save
# plt.tight_layout()
# plt.savefig('pi_visualization.png')

# # Print summary statistics
# print("\nSummary Statistics:")
# print(df[['log_pi_present', 'log_pi_absent', 'log_pi_between']].describe())