#!/usr/bin/env python3
import argparse
import gffutils
import logging
from collections import defaultdict

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_db(gff_file, is_miniprot=True):
    """Create a GFF database"""
    
    # Transform function to handle miniprot's attributes
    def transform_func(x):
        # Keep all attributes
        return x

    try:
        db = gffutils.create_db(
            gff_file,
            dbfn=':memory:',
            force=True,
            merge_strategy='create_unique',
            transform=transform_func,
            disable_infer_genes=True,
            disable_infer_transcripts=True
        )
        return db
    except Exception as e:
        logger.error(f"Error creating database: {str(e)}")
        raise

def extract_ref_gene_info(ref_db):
    """Extract reference gene information"""
    ref_genes = {}
    for gene in ref_db.features_of_type('gene'):
        gene_id = gene.id
        ref_genes[gene_id] = {
            'length': gene.end - gene.start + 1,
            'coords': (gene.start, gene.end)
        }
    return ref_genes

def get_valid_features(db, min_identity=0.7):
    """Get valid features based on identity threshold"""
    valid_features = set()
    kept_mrnas = set()
    
    # First process mRNA features and their identity scores
    for feature in db.features_of_type('mRNA'):
        try:
            identity = float(feature.attributes.get('Identity', [0])[0])
            if identity >= min_identity:
                kept_mrnas.add(feature.id)
                valid_features.add(feature)
                # Also keep all child features (CDS)
                for child in db.children(feature):
                    valid_features.add(child)
        except ValueError as e:
            logger.warning(f"Could not parse identity score for feature {feature.id}: {e}")
            continue
    
    return valid_features

def filter_gff(input_gff, output_gff, min_identity=0.7):
    """Filter GFF based on identity threshold"""
    logger.info("Creating database...")
    db = create_db(input_gff, is_miniprot=True)
    
    logger.info("Finding valid features...")
    valid_features = get_valid_features(db, min_identity)
    
    logger.info(f"Found {len(valid_features)} valid features")
    
    # Write filtered GFF
    logger.info("Writing filtered GFF...")
    with open(output_gff, 'w') as out:
        # Write GFF version header
        out.write("##gff-version 3\n")
        
        # Write valid features
        for feature in valid_features:
            out.write(str(feature) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Filter miniprot GFF annotations')
    parser.add_argument('--input', required=True, help='Input GFF from miniprot')
    parser.add_argument('--output', required=True, help='Output filtered GFF')
    parser.add_argument('--min-identity', type=float, default=0.7,
                       help='Minimum protein identity for feature validation')
    
    args = parser.parse_args()
    
    try:
        filter_gff(args.input, args.output, args.min_identity)
    except Exception as e:
        logger.error(f"Error during filtering: {str(e)}")
        raise

if __name__ == '__main__':
    main()