import pandas as pd
import gffutils
from pathlib import Path
import argparse

def load_and_filter_calls(calls_file):
    """Load calls matrix and identify polymorphic ortholog groups."""
    # Read calls matrix
    calls_df = pd.read_csv(calls_file, sep='\t')
    
    # Drop RCC1749 and RCC3052 columns
    cols_to_drop = ['RCC1749', 'RCC3052']
    filtered_df = calls_df.drop(columns=cols_to_drop)
    
    # Function to check if row is polymorphic (contains both 0 and 1, no 3s)
    def is_polymorphic(row):
        values = row.drop('ortholog_id')
        return (values.isin([0, 1]).all() and  # only contains 0s and 1s
                (values == 0).any() and        # has at least one 0
                (values == 1).any())           # has at least one 1
    
    # Filter for polymorphic rows
    polymorphic_df = filtered_df[filtered_df.apply(is_polymorphic, axis=1)]
    
    return polymorphic_df['ortholog_id'].tolist()

def load_metadata(metadata_file, polymorphic_ids):
    """Load and filter metadata for polymorphic ortholog groups."""
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    # Filter for polymorphic IDs and exclude RCC1749 and RCC3052
    metadata_df = metadata_df[
        (metadata_df['ortholog_id'].isin(polymorphic_ids)) & 
        (~metadata_df['sample_name'].isin(['RCC1749', 'RCC3052']))
    ]
    return metadata_df

def adjust_coordinates(row):
    """Adjust coordinates by trimming 100bp flanking regions."""
    return pd.Series({
        'adj_start': int(row['start']) + 100 if row['start'] != '?' else None,
        'adj_end': int(row['end']) - 100 if row['end'] != '?' else None
    })

def find_gene_associations(gtf_db, contig, start, end, upstream_dist=1000):
    """Find overlapping or nearby genes for a given genomic region."""
    # First check for direct overlaps
    overlapping = list(gtf_db.region(seqid=contig, start=start, end=end, featuretype='gene'))
    if overlapping:
        return [(gene['gene_id'][0], 'overlapping', 0) for gene in overlapping]
    
    # If no overlaps, check upstream/downstream
    nearby_genes = []
    
    # Check upstream
    upstream_genes = list(gtf_db.region(
        seqid=contig,
        start=max(0, start - upstream_dist),
        end=start,
        featuretype='gene'
    ))
    for gene in upstream_genes:
        distance = start - gene.end
        if distance <= upstream_dist:
            nearby_genes.append((gene['gene_id'][0], 'upstream', distance))
    
    # Check downstream
    downstream_genes = list(gtf_db.region(
        seqid=contig,
        start=end,
        end=end + upstream_dist,
        featuretype='gene'
    ))
    for gene in downstream_genes:
        distance = gene.start - end
        if distance <= upstream_dist:
            nearby_genes.append((gene['gene_id'][0], 'downstream', distance))
    
    return nearby_genes if nearby_genes else [('No_gene', 'none', None)]

def create_or_get_gtf_db(sample_name, gtf_dir):
    """Create or get GTF database for a specific sample."""
    gtf_file = Path(gtf_dir) / f"{sample_name}.vg_paths.gtf"
    db_file = Path(gtf_dir) / f"{sample_name}.vg_paths.gtf.db"
    
    if not db_file.exists():
        print(f"Creating GTF database for {sample_name}...")
        db = gffutils.create_db(
            str(gtf_file),
            str(db_file),
            merge_strategy='merge',
            force=True,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
    return gffutils.FeatureDB(str(db_file), keep_order=True)

def main(calls_file, metadata_file, gtf_dir, output_file):
    # Dictionary to store GTF databases for each sample
    gtf_dbs = {}
    
    # Get polymorphic ortholog IDs
    polymorphic_ids = load_and_filter_calls(calls_file)
    print(f"Found {len(polymorphic_ids)} polymorphic ortholog groups")
    
    # Load and filter metadata
    metadata_df = load_metadata(metadata_file, polymorphic_ids)
    
    # Adjust coordinates
    metadata_df[['adj_start', 'adj_end']] = metadata_df.apply(adjust_coordinates, axis=1)
    
    # Initialize results list
    results = []
    
    # Process each polymorphic ortholog
    for _, row in metadata_df.iterrows():
        if row['contig'] != '?' and row['adj_start'] is not None:
            # Get or create sample-specific GTF database
            sample_name = row['sample_name']
            if sample_name not in gtf_dbs:
                gtf_dbs[sample_name] = create_or_get_gtf_db(sample_name, gtf_dir)
            
            gene_associations = find_gene_associations(
                gtf_dbs[sample_name],
                row['contig'],
                row['adj_start'],
                row['adj_end']
            )
            
            for gene_id, relationship, distance in gene_associations:
                results.append({
                    'ortholog_id': row['ortholog_id'],
                    'sample': row['sample_name'],
                    'contig': row['contig'],
                    'start': row['start'],
                    'end': row['end'],
                    'gene_id': gene_id,
                    'relationship': relationship,
                    'distance': distance
                })
    
    # Convert results to DataFrame and save
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Analyze polymorphic introner elements')
    parser.add_argument('--calls', required=True, help='Path to calls matrix file')
    parser.add_argument('--metadata', required=True, help='Path to metadata file')
    parser.add_argument('--gtf-dir', required=True, help='Directory containing sample GTF files')
    parser.add_argument('--output', required=True, help='Path to output file')
    
    args = parser.parse_args()
    
    main(args.calls, args.metadata, args.gtf_dir, args.output)