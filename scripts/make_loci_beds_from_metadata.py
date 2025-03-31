import pandas as pd
import sys

sample_name = sys.argv[2]

tsv_file = sys.argv[1]

output_file = sys.argv[3]

# Read the TSV file into a DataFrame
df = pd.read_csv(tsv_file, sep='\t', header=0)

# Filter rows for the specified sample_name
sample_df = df[df['sample_name'] == sample_name]

# Check if any data exists for the sample
if sample_df.empty:
    print(f"No data found for sample: {sample_name}")
    sys.exit(1)

sample_df = sample_df.replace('?', pd.NA).dropna()
sample_df[["start", "end"]] = sample_df[["start", "end"]].astype(int)

# Sort by 'contig' and 'start'
sample_sorted = sample_df.sort_values(['contig', 'start'])

# Select BED-specific columns: 'contig', 'start', 'end'
bed_df = sample_sorted[['contig', 'start', 'end', 'sample_name', 'ortholog_id']]

# Define the output filename
output_filename = output_file

# Write to the BED file without headers or index
bed_df.to_csv(output_filename, sep='\t', header=False, index=False)

print(f"BED file created: {output_filename}")