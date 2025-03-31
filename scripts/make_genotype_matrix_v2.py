import pandas as pd
import sys

def split_tsv(input_tsv, output_file1, output_file2):
    # Load TSV file
    df = pd.read_csv(input_tsv, sep="\t")

    # **File 1: Keep necessary columns**
    df_file1 = df[['ortholog_id', 'sample', 'contig', 'left_flank_start', 'right_flank_end']]
    df_file1.rename(columns={'sample':'sample_name', 'left_flank_start':'start', 'right_flank_end':'end'}, inplace=True)
    df_file1.to_csv(output_file1, sep="\t", index=True)

    # **File 2: Pivot to wide format (ortholog_id as rows, samples as columns)**
    df_file2 = df.pivot(index='ortholog_id', columns='sample', values='presence')

    # Reset index to make sure `ortholog_id` remains a column
    df_file2.reset_index(inplace=True)

    # Save to file
    df_file2.to_csv(output_file2, sep="\t", index=False)

    print(f"Metadata file saved: {output_file1}")
    print(f"Calls file saved: {output_file2}")

if __name__ == "__main__":
    """
    Command-line usage:
    python split_tsv.py <input_tsv> <output_file1> <output_file2>

    Example:
    python split_tsv.py my_data.tsv file1.tsv file2.tsv
    """
    if len(sys.argv) != 4:
        print("Usage: python split_tsv.py <input_tsv> <output_file1> <output_file2>")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_file1 = sys.argv[2]
    output_file2 = sys.argv[3]

    split_tsv(input_tsv, output_file1, output_file2)