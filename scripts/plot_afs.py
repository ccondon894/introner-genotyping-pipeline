import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LogNorm
call_in = sys.argv[1]

figure_afs_2d = sys.argv[2]
figure_afs_1 = sys.argv[3]
figure_afs_2 = sys.argv[4]

df = pd.read_csv(call_in, sep="\t", header=0, index_col=0)

group2 = ['RCC1749', 'RCC3052']
group1 = [col for col in df.columns if col not in group2]

n_group1= len(group1)
n_group2 = len(group2)


# Initialize the 2D AFS matrix
afs_2d = np.zeros((n_group1 + 1, n_group2 + 1), dtype=int)

# Initialize 1D AFS arrays for each group
afs1 = np.zeros(n_group1 + 1, dtype=int)
afs2 = np.zeros(n_group2 + 1, dtype=int)

# Step 3: Iterate over each locus/row in the DataFrame
for idx, row in df.iterrows():
    # Extract genotype calls for each group at the current locus
    group1_calls = row[group1]
    group2_calls = row[group2]
    
    # ---------- 2D AFS Update ----------
    # Check for missing data in either group for joint analysis
    if not ((group1_calls == 3).any() or (group2_calls == 3).any()):
        # Count "present" calls for joint analysis
        count1 = (group1_calls == 0).sum()
        count2 = (group2_calls == 0).sum()
        afs_2d[count1, count2] += 1

    # ---------- 1D AFS for Group 1 ----------
    # Check for missing data in group1 only
    if not (group1_calls == 3).any():
        # Count "present" calls in group1
        count1_only = (group1_calls == 0).sum()
        afs1[count1_only] += 1

    # ---------- 1D AFS for Group 2 ----------
    # Check for missing data in group2 only
    if not (group2_calls == 3).any():
        # Count "present" calls in group2
        count2_only = (group2_calls == 0).sum()
        afs2[count2_only] += 1

# Now, afs1 and afs2 include bins from 0 up to n_group
# We can simply ignore or remove the 0-count bin for reporting

# Exclude the first bin (index 0) for reporting:
afs1_no_zero = afs1[1:-1]  # This slices out the first element
afs2_no_zero = afs2[1:]

# print("1D AFS for Group 1 (ignoring loci with 0 present calls):")
# for allele_count, frequency in enumerate(afs1_no_zero, start=1):
#     print(f"Allele count {allele_count}: {frequency}")

# print("\n1D AFS for Group 2 (ignoring loci with 0 present calls):")
# for allele_count, frequency in enumerate(afs2_no_zero, start=1):
#     print(f"Allele count {allele_count}: {frequency}")


print(afs_2d)

print(afs1_no_zero)

print(afs2_no_zero)

log10_afs_2d = np.log10(afs_2d + 1)
log10_afs_2d = log10_afs_2d.T

plt.figure(figsize=(20, 5))

sns.heatmap(log10_afs_2d, fmt=".1f", annot=True, cmap="viridis", cbar_kws={'label': 'log10 Introner Count', "shrink": 1.0}, square=True, linewidths=0.5, linecolor='white')
plt.ylabel('Allele Count in Group 2 (size 2)')
plt.xlabel('Allele Count in Group 1 (size 14)')
plt.gca().invert_yaxis()
plt.tight_layout()
plt.savefig(figure_afs_2d)

plt.figure(figsize=(15,10))
plt.bar(range(1, len(afs1_no_zero)+1), afs1_no_zero, color='skyblue')
xticks = range(0, len(afs1_no_zero)+1)
plt.xticks(ticks=np.arange(len(xticks)), labels=xticks, fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Frequency', fontsize=14)
plt.xlabel('Allele Count Group 1', fontsize=14)
plt.tight_layout()
plt.savefig(figure_afs_1)

plt.figure(figsize=(5,10))
plt.bar(range(1, len(afs2_no_zero) + 1), afs2_no_zero, color='skyblue')
xticks = range(0, len(afs2_no_zero) + 1)
plt.xticks(ticks=np.arange(len(xticks)), labels=xticks, fontsize=14)
plt.yticks(fontsize=14)
plt.ylabel('Frequency',fontsize=14)
plt.xlabel('Allele Count Group 2', fontsize=14)
plt.tight_layout()
plt.savefig(figure_afs_2)