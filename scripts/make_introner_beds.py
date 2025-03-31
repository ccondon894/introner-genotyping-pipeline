import sys
from collections import defaultdict
import pandas as pd

md_file = sys.argv[1]
calls_file = sys.argv[2]

introner_dict = defaultdict(list)
calls_df = pd.read_csv(calls_file, sep="\t")

with open(md_file, 'r') as f:
    next(f)
    for line in f:
        line = line.strip().split("\t")
        line = line[1:]
        sample = line[1]
        ortho_id = line[0]
        # print(sample, ortho_id)
        # print(calls_df.loc[(calls_df['ortholog_id'] == ortho_id), sample])
        if line[2] != '?' and int(calls_df.loc[(calls_df['ortholog_id'] == ortho_id), sample]) != 1:
            contig = line[2]
            start = int(line[3]) + 100
            end = int(line[4]) - 100

            region = [contig, start, end]
            introner_dict[sample].append(region)

introner_dict = dict(introner_dict)
for sample in introner_dict:
    outfile = f"/scratch1/chris/introner-genotyping-pipeline/introner_beds/{sample}.introner_loci.bed"
    regions = introner_dict[sample]
    sorted_regions = sorted(regions, key=lambda x: (x[0], int(x[1]), int(x[2])))
    with open(outfile, 'w') as o:
        for region in sorted_regions:
            region_str = f"{region[0]}\t{region[1]}\t{region[2]}\n"
            o.write(region_str)
        
        o.close()
    print(f"Wrote {sample} loci to {outfile}")