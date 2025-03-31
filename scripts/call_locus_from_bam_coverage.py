import pandas as pd
import sys

bed_file = sys.argv[1]
depth_file = sys.argv[2]
call_file = sys.argv[3]
output_bed_file = sys.argv[4]

depth_df = pd.read_csv(depth_file, sep="\t", header=None, names=["chrom", "pos", "depth"])
depth_df =depth_df.sort_values(["chrom", "pos"]).reset_index(drop=True)

call_df = pd.read_csv(call_file, sep="\t", header=0, index_col=0)

bed_df = pd.read_csv(bed_file, sep="\t", header=None, names=["chrom", "start", "end", "sample_id", "ortholog_id"])
# Function to determine presence/absence for a single interval
def call_presence(row, call, depth_df):
    chrom = row["chrom"]
    start = int(row["start"])
    end = int(row["end"])
    
    left_start, left_end = start, start + 100
    right_start, right_end = end - 100, end
    # Load coverage for this interval from depth_df
    # Filter depth_df for this chromosome and region
    interval_cov = depth_df[(depth_df["chrom"] == chrom) & (depth_df["pos"] >= start) & (depth_df["pos"] < end)]
    # Check the left and right regions
    left_cov = interval_cov[(interval_cov["pos"] >= left_start) & (interval_cov["pos"] < left_end)]
    right_cov = interval_cov[(interval_cov["pos"] >= right_start) & (interval_cov["pos"] < right_end)]

    # Left region
    left_sum_depth = left_cov["depth"].sum()
    # missing bases in left region = zero depth
    left_mean_depth = left_sum_depth / 100  # Always 100 bases expected
    # Right region
    right_sum_depth = right_cov["depth"].sum()
    # missing bases in right region = zero depth
    right_mean_depth = right_sum_depth / 100

    if call == 1:
        new_call = None
        mean_middle_depth = 0
        # Presence/Absence logic
        if left_mean_depth < 10 or right_mean_depth < 10: # insufficient coverage on left/right flanks
            new_call = 3
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call
        
        if not new_call:
            return left_mean_depth, mean_middle_depth, right_mean_depth, call

    else:
        middle_start, middle_end = left_end, right_start

        # Junction points
        junction_LM = middle_start   # boundary between left and middle
        junction_MR = middle_end     # boundary between middle and right
        # Check no zero-depth regions near junctions:
        # For L-M junction: range = [junction_LM-10, junction_LM+10)
        # For M-R junction: range = [junction_MR-10, junction_MR+10)
        LM_range = interval_cov[(interval_cov["pos"] >= junction_LM - 10) & (interval_cov["pos"] < junction_LM + 10)]
        MR_range = interval_cov[(interval_cov["pos"] >= junction_MR - 10) & (interval_cov["pos"] < junction_MR + 10)]

        # Check middle region criteria:
        middle_cov = interval_cov[(interval_cov["pos"] >= middle_start) & (interval_cov["pos"] < middle_end)]
        middle_len = middle_end - middle_start
        # If we don't have complete coverage rows for the entire middle region, assume missing positions have zero depth
        # Identify which positions in the range [middle_start, middle_end) are missing
        covered_positions = set(middle_cov["pos"].unique())
        all_middle_positions = set(range(middle_start, middle_end))
        missing_positions = all_middle_positions - covered_positions
        # Create a "complete" coverage array for the middle region
        # Positions not listed in middle_cov are implicitly zero depth

        # For performance, we'll just measure counts and sums:
        nonzero_count = (middle_cov["depth"] > 0).sum()
        total_nonzero = nonzero_count  # from covered positions

        # no coverage info means depth=0, so no increment
        total_covered_bases = middle_len
        total_nonzero += 0  # no coverage means zero depth

        sum_middle_depth = middle_cov["depth"].sum()
        # Add zeros for missing positions
        sum_middle_depth += 0 * len(missing_positions)  # no addition needed, just clarifies logic
        mean_middle_depth = sum_middle_depth / total_covered_bases

        new_call = None
        # Presence/Absence logic
        if left_mean_depth < 10 or right_mean_depth < 10: # insufficient coverage on left/right flanks
            new_call = 3
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call

        if (mean_middle_depth > 5 * left_mean_depth) or (mean_middle_depth > 5 * right_mean_depth): # insufficient coverage on left/right flanks
            new_call = 3
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call

        # At least 80% of the bases in the middle have non-zero depth
        if total_nonzero < 0.8 * total_covered_bases:
            new_call = 1
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call
        
        # Mean depth of the middle region >= 10
        if mean_middle_depth < 10:
            new_call = 1
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call

        # The mean depth of the left or right regions must not exceed 5X the mean depth of the middle region
        if (left_mean_depth > 5 * mean_middle_depth) or (right_mean_depth > 5 * mean_middle_depth):
            new_call = 1
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call

        if LM_range.shape[0] < 20 or MR_range.shape[0] < 20:
            # This means we don't have coverage info for all those bases (some missing?), 
            # if missing positions exist, consider them zero coverage. So fail if not fully covered:
            new_call = 1
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call

        if (LM_range["depth"] == 0).any() or (MR_range["depth"] == 0).any():
            new_call = 1
            return left_mean_depth, mean_middle_depth, right_mean_depth, new_call

        if not new_call: # if all criteria are met, then we just return the original call 
            return left_mean_depth, mean_middle_depth, right_mean_depth, call

with open(output_bed_file, 'w') as o:
    o.write(f"{'chrom'}\t{'start'}\t{'end'}\t{'sample_id'}\t{'ortholog_id'}\t{'left_mean_depth'}\t{'mean_middle_depth'}\t{'right_mean_depth'}\t{'call'}\t{'new_call'}\n")
    for index, row in bed_df.iterrows():
        print("Checking locus...", row["ortholog_id"], row["sample_id"])
        call = call_df.loc[row["ortholog_id"], row["sample_id"]]
        if call == 3:
            print("skipping locus due to missing data")

        else:
            left_mean_depth, mean_middle_depth, right_mean_depth, new_call = call_presence(row, call, depth_df)
            data_row = f"{row['chrom']}\t{row['start']}\t{row['end']}\t{row['sample_id']}\t{row['ortholog_id']}\t{left_mean_depth}\t{mean_middle_depth}\t{right_mean_depth}\t{call}\t{new_call}\n"
            o.write(data_row)
            print(left_mean_depth, mean_middle_depth, right_mean_depth, call, new_call, sep="\t")
            
