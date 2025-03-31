# combine_ortholog_results.py
import argparse
import pandas as pd

def main(args):
    dfs = []
    for file in args.input_files:
        sample = file.split('/')[-1].split('.')[0].split('_')[1]
        df = pd.read_csv(file, sep='\t')
        df['sample'] = sample
        dfs.append(df)
    
    pd.concat(dfs).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_files', nargs='+')
    parser.add_argument('--output')
    main(parser.parse_args())