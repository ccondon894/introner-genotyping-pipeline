import glob
import sys
import os
import pandas as pd
from Bio import SeqIO

# Each ortholog group will contain the following files:
# {ortholog_id}.{left}_flank.present.fa
# {ortholog_id}.{right}_flank.present.fa
# {ortholog_id}.{left}_flank.absent.fa
# {ortholog_id}.{right}_flank.absent.fa

metadata_file = sys.argv[1]
calls_file = sys.argv[2]
fasta_dir = sys.argv[3]
output_dir = sys.argv[4]

meta_df = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)
meta_g1 = meta_df[~meta_df['sample_name'].isin(['RCC1749', 'RCC3052'])]

calls_df = pd.read_csv(calls_file, sep="\t", header=0)
calls_g1 = calls_df.drop(['RCC1749', 'RCC3052'], axis=1)
calls_g1 = calls_g1[~calls_g1.isin([3]).any(axis=1)]

fasta_dict = {}
for file in glob.glob(f"{fasta_dir}/*.consensus.fa"):
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_dict[record.id] = str(record.seq)


os.makedirs(output_dir, exist_ok=True)

for _, row in calls_g1.iterrows():

    ortholog_id = row['ortholog_id']
    present_samples = [col for col in row.index if row[col] == 0]
    absent_samples = [col for col in row.index if row[col] != 0 and col != 'ortholog_id']
    # print(present_samples)
    # print(absent_samples)
    # skip any ortholog groups where 
    if len(present_samples) < 2 or len(absent_samples) < 2:
        continue

    present_left_file = f"{output_dir}/{ortholog_id}.left_flank.present.fa"
    present_right_file = f"{output_dir}/{ortholog_id}.right_flank.present.fa"

    absent_left_file = f"{output_dir}/{ortholog_id}.left_flank.absent.fa"
    absent_right_file = f"{output_dir}/{ortholog_id}.right_flank.absent.fa"

    with open(present_left_file, 'w') as l, open(present_right_file, 'w') as r:

        for sample in present_samples:
            coords_row = meta_g1[(meta_g1['ortholog_id'] == ortholog_id) & (meta_g1['sample_name'] == sample)]
            left_start, left_end = int(coords_row['start'].iloc[0]), int(coords_row['start'].iloc[0]) + 100
            right_start, right_end = int(coords_row['end'].iloc[0]) - 100, int(coords_row['end'].iloc[0])

            left_flank_id = f"{coords_row['contig'].iloc[0]}:{left_start}-{left_end}"
            left_flank_seq = fasta_dict[left_flank_id]
            
            right_flank_id = f"{coords_row['contig'].iloc[0]}:{right_start}-{right_end}"
            right_flank_seq = fasta_dict[right_flank_id]

            l.write(f">{left_flank_id}\n")
            l.write(f"{left_flank_seq}\n")

            r.write(f">{right_flank_id}\n")
            r.write(f"{right_flank_seq}\n")
        

    with open(absent_left_file, 'w') as l, open(absent_right_file, 'w') as r:

        for sample in absent_samples:
            coords_row = meta_g1[(meta_g1['ortholog_id'] == ortholog_id) & (meta_g1['sample_name'] == sample)]
            left_start, left_end = int(coords_row['start'].iloc[0]), int(coords_row['start'].iloc[0]) + 100
            right_start, right_end = int(coords_row['end'].iloc[0]) - 100, int(coords_row['end'].iloc[0])

            left_flank_id = f"{coords_row['contig'].iloc[0]}:{left_start}-{left_end}"
            left_flank_seq = fasta_dict[left_flank_id]
            
            right_flank_id = f"{coords_row['contig'].iloc[0]}:{right_start}-{right_end}"
            right_flank_seq = fasta_dict[right_flank_id]

            l.write(f">{left_flank_id}\n")
            l.write(f"{left_flank_seq}\n")

            r.write(f">{right_flank_id}\n")
            r.write(f"{right_flank_seq}\n")

        


