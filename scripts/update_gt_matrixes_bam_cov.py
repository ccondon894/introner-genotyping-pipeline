import pandas as pd
import sys
import os


metadata_in = sys.argv[1]
call_in = sys.argv[2]
bed_path = sys.argv[3]

metadata_out = sys.argv[4]
call_out = sys.argv[5]

call_df = pd.read_csv(call_in, sep="\t", header=0, index_col=0)
meta_df = pd.read_csv(metadata_in, sep="\t", header=0, index_col=0)

# update genotype matrix with new calls
for file in os.listdir(bed_path):
    if file.endswith("coverage_calls.bed"):
        print("Reading...", file)
        file_path = os.path.join(bed_path, file)
        sample = file.split(".")[0]
        bed_df = pd.read_csv(file_path, sep="\t", header=0)

        for index, row in bed_df.iterrows():
            chrom, start, end, ortho_id, new_call = row['chrom'], row['start'], row['end'], row['ortholog_id'], int(row['new_call'])
            # update call df with new call
            call_df.loc[ortho_id, sample] = new_call


for index, row in meta_df.iterrows():
    sample_id = row['sample_name']
    ortho_id = row['ortholog_id']
    if call_df.loc[ortho_id, sample_id] == 3:
        meta_df.at[index, 'contig'] = '?'
        meta_df.at[index, 'start'] = '?'
        meta_df.at[index, 'end'] = '?'


call_df.to_csv(call_out, sep="\t", columns=call_df.columns)
meta_df.to_csv(metadata_out, sep="\t", columns=meta_df.columns)
