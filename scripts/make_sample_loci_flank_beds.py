import sys
import os
import pandas as pd


metadata_file = sys.argv[1]
sample_name = sys.argv[2]
left_bed = sys.argv[3]
right_bed = sys.argv[4]


meta_df = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)

sample_df = meta_df[meta_df['sample_name'] == sample_name]
sample_df = sample_df[sample_df['start'] != '?']
sample_df['start'] = sample_df['start'].astype(int)
sample_df['end'] = sample_df['end'].astype(int)
sample_df = sample_df.sort_values(by=['contig', 'start'])

with open(left_bed, 'w') as l, open(right_bed, 'w') as r:
    for index, row in sample_df.iterrows():

        contig = row['contig']
        left_start = row['start']
        left_end = row['start'] + 100
        right_start = row['end'] - 100
        right_end = row['end']
        ortholog_id = row['ortholog_id']

        l.write(f"{contig}\t{left_start}\t{left_end}\t{ortholog_id}\n")
        r.write(f"{contig}\t{right_start}\t{right_end}\t{ortholog_id}\n")