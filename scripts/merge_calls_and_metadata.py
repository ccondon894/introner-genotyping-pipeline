import pandas as pd
import sys
# Load the files
calls_path = sys.argv[1]
metadata_path = sys.argv[2]
output_path = sys.argv[3]

# Read the files
calls_df = pd.read_csv(calls_path, sep="\t")
metadata_df = pd.read_csv(metadata_path, sep="\t")

# Create a dictionary to store calls for each ortholog_id and sample
call_dict = {}
samples = calls_df.columns[1:]  # Skip the ortholog_id column

# Populate the dictionary with calls
for _, row in calls_df.iterrows():
    ortholog_id = row['ortholog_id']
    for sample in samples:
        key = (ortholog_id, sample)
        call_dict[key] = row[sample]

# Add a 'call' column to the metadata DataFrame
metadata_df['call'] = metadata_df.apply(
    lambda row: call_dict.get((row['ortholog_id'], row['sample_name']), None), 
    axis=1
)

# Write to TSV
metadata_df.to_csv(output_path, sep="\t", index=False)

print(f"Merged matrix saved to {output_path}")