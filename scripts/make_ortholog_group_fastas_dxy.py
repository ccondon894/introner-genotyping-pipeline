import glob
import sys
import os
import pandas as pd
from Bio import SeqIO

# Each ortholog group will contain the following files:
# {ortholog_id}.g1.{left}_flank.{presence}.fa
# {ortholog_id}.g1.{right}_flank.{presence}.fa
# {ortholog_id}.g2.{left}_flank.{presence}.fa
# {ortholog_id}.g2.{right}_flank.{presence}.fa

metadata_file = sys.argv[1]
calls_file = sys.argv[2]
fasta_dir = sys.argv[3]
output_dir = sys.argv[4]

meta_df = pd.read_csv(metadata_file, sep="\t", header=0, index_col=0)

calls_df = pd.read_csv(calls_file, sep="\t", header=0)
calls_df = calls_df[~calls_df.isin([3]).any(axis=1)]

fasta_dict = {}
for file in glob.glob(f"{fasta_dir}/*.consensus.fa"):
    with open(file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_dict[record.id] = str(record.seq)


os.makedirs(output_dir, exist_ok=True)

for _, row in calls_df.iterrows():

    ortholog_id = row['ortholog_id']
    group1  = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC2482", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]
    group2 = ["RCC1749", "RCC3052"]
    present_group1 = [col for col in group1 if row[col] == 0]
    present_group2 = [col for col in group2 if row[col] == 0]
    presence_g1 = None
    presence_g2 = None

    if len(present_group1) == 14 and len(present_group2) == 2:
        presence_g1 = "present"
        presence_g2 = "present"
    
    elif len(present_group1) == 14 and len(present_group2) == 0:
        presence_g1 = "present"
        presence_g2 = "absent"
    
    elif len(present_group1) == 0 and len(present_group2) == 2:
        presence_g1 = "absent"
        presence_g2 = "present"

    else:
        continue


    g1_left_file = f"{output_dir}/{ortholog_id}.g1.left_flank.{presence_g1}.fa"
    g1_right_file = f"{output_dir}/{ortholog_id}.g1.right_flank.{presence_g1}.fa"

    g2_left_file = f"{output_dir}/{ortholog_id}.g2.left_flank.{presence_g2}.fa"
    g2_right_file = f"{output_dir}/{ortholog_id}.g2.right_flank.{presence_g2}.fa"

    with open(g1_left_file, 'w') as l, open(g1_right_file, 'w') as r:

        for sample in group1:
            coords_row = meta_df[(meta_df['ortholog_id'] == ortholog_id) & (meta_df['sample_name'] == sample)]
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
        

    with open(g2_left_file, 'w') as l, open(g2_right_file, 'w') as r:

        for sample in group2:
            coords_row = meta_df[(meta_df['ortholog_id'] == ortholog_id) & (meta_df['sample_name'] == sample)]
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

        


