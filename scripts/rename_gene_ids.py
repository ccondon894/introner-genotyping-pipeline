#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GTF attribute string into a dictionary"""
    attrs = {}
    # Split by semicolon and strip whitespace
    for attr in attr_str.strip().split(';'):
        if attr.strip():
            try:
                # Split by first space and strip quotes
                key, value = attr.strip().split(' ', 1)
                attrs[key] = value.strip('"')
            except ValueError:
                continue
    return attrs

def format_attributes(attrs):
    """Format attribute dictionary back to GTF format"""
    formatted = []
    # Ensure gene_id and transcript_id come first if they exist
    priority_attrs = ['gene_id', 'transcript_id']
    for key in priority_attrs:
        if key in attrs:
            formatted.append(f'{key} "{attrs[key]}"')
    
    # Add remaining attributes
    for key, value in attrs.items():
        if key not in priority_attrs:
            formatted.append(f'{key} "{value}"')
    
    return '; '.join(formatted) + ';'

def extract_target_id(attrs):
    """Extract target ID from attributes"""
    if 'Target' in attrs:
        return attrs['Target'].split()[0]
    return None

def create_reference_mappings(ref_gtf):
    """Create mappings from target IDs to reference gene/transcript IDs"""
    transcript_map = {}
    gene_map = {}
    
    with open(ref_gtf) as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) != 9:
                continue
                
            attrs = parse_attributes(fields[8])
            
            # For transcript lines, store the mapping
            if fields[2] == 'transcript':
                if 'transcript_id' in attrs:
                    transcript_id = attrs['transcript_id']
                    gene_id = attrs.get('gene_id')
                    transcript_map[transcript_id] = {
                        'transcript_id': transcript_id,
                        'gene_id': gene_id
                    }
    
    return transcript_map

def update_gtf(input_gtf, output_gtf, transcript_map):
    """Update query GTF with reference IDs"""
    with open(input_gtf) as f, open(output_gtf, 'w') as out:
        for line in f:
            if line.startswith('#'):
                out.write(line)
                continue
                
            fields = line.strip().split('\t')
            if len(fields) != 9:
                out.write(line)
                continue
            
            attrs = parse_attributes(fields[8])
            target_id = extract_target_id(attrs)
            
            if target_id and target_id in transcript_map:
                ref_ids = transcript_map[target_id]
                
                # Update gene_id and transcript_id if present
                if 'gene_id' in attrs:
                    attrs['gene_id'] = ref_ids['gene_id']
                if 'transcript_id' in attrs and fields[2] != 'gene':
                    attrs['transcript_id'] = ref_ids['transcript_id']
                if 'Parent' in attrs:
                    if fields[2] in ['mRNA', 'transcript']:
                        attrs['Parent'] = ref_ids['gene_id']
                    else:
                        attrs['Parent'] = ref_ids['transcript_id']
                
                # Update the ID field if present
                if 'ID' in attrs:
                    if fields[2] == 'gene':
                        attrs['ID'] = ref_ids['gene_id']
                    elif fields[2] in ['mRNA', 'transcript']:
                        attrs['ID'] = ref_ids['transcript_id']
                    # Keep other IDs as they are (exons, CDS)
            
            # Update the attributes field
            fields[8] = format_attributes(attrs)
            out.write('\t'.join(fields) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Update query GTF IDs to match reference')
    parser.add_argument('--reference', required=True, help='Reference GTF file')
    parser.add_argument('--query', required=True, help='Query GTF file to update')
    parser.add_argument('--output', required=True, help='Output GTF file')
    
    args = parser.parse_args()
    
    # Create mappings from reference GTF
    print("Creating reference mappings...")
    transcript_map = create_reference_mappings(args.reference)
    print(f"Found {len(transcript_map)} reference transcripts")
    
    # Update query GTF
    print("Updating query GTF...")
    update_gtf(args.query, args.output, transcript_map)
    print("Done!")

if __name__ == '__main__':
    main()