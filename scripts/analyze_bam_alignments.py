import argparse
import logging
import pysam
from Bio import SeqIO
import pandas as pd

def parse_best_alignments(bam_file):
    """Extract best alignments for each query sequence based on MAPQ and CIGAR matches."""
    alignments = {}
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    for read in bam:
        query_name = read.query_name
        
        # Skip unmapped reads
        if read.is_unmapped:
            continue
            
        # Count matches in CIGAR string
        matches = sum([length for op, length in read.cigartuples if op == 0])
        
        # Store alignment if it's better than current best
        if query_name not in alignments or \
           (read.mapping_quality > alignments[query_name]['mapq']) or \
           (read.mapping_quality == alignments[query_name]['mapq'] and 
            matches > alignments[query_name]['matches']):
            alignments[query_name] = {
                'contig': bam.get_reference_name(read.reference_id),
                'start': read.reference_start,
                'end': read.reference_end,
                'mapq': read.mapping_quality,
                'matches': matches,
                'is_reverse': read.is_reverse
            }
    
    bam.close()
    return alignments

def find_ortholog_pairs(left_alignments, right_alignments):
    """Find valid left-right pairs and determine ortholog coordinates."""
    ortholog_pairs = []
    
    for query_name in left_alignments:
        base_id = query_name.replace('_left', '')
        right_name = base_id + '_right'
        
        left_aln = left_alignments[query_name]
        
        # Case: missing right flank alignment
        if right_name not in right_alignments:
            ortholog_pairs.append((base_id, left_aln['contig'], 
                                 left_aln['start'], None, 3))
            continue
        
        right_aln = right_alignments[right_name]
        
        # Check if alignments are on same contig
        if left_aln['contig'] != right_aln['contig']:
            ortholog_pairs.append((base_id, left_aln['contig'], 
                                 left_aln['start'], None, 3))
            continue
            
        # Get coordinates considering strand
        left_end = max(left_aln['start'], left_aln['end'])
        right_start = min(right_aln['start'], right_aln['end'])
        
        # Check distance between flanks
        distance = right_start - left_end
        if distance > 300 or right_aln['end'] <= left_end:
            ortholog_pairs.append((base_id, left_aln['contig'], 
                                 left_aln['start'], None, 3))
            continue
            
        ortholog_pairs.append((base_id, left_aln['contig'],
                             min(left_aln['start'], left_aln['end']),
                             max(right_aln['start'], right_aln['end']),
                             None))
    
    return ortholog_pairs

def check_locus_presence(query_fa, target_fa, ortholog_pairs):
    """Check for presence/absence of introner at each locus."""
    query_seqs = {record.id: str(record.seq) for record in SeqIO.parse(query_fa, 'fasta')}
    target_seqs = {record.id: str(record.seq) for record in SeqIO.parse(target_fa, 'fasta')}
    
    results = []
    for query_id, target_contig, left_pos, right_pos, presence in ortholog_pairs:
        if presence == 3:  # Missing flanks
            results.append([query_id, target_contig, left_pos, right_pos, 3])
            continue
            
        query_seq = query_seqs[query_id]
        target_seq = target_seqs[target_contig][left_pos:right_pos]
        
        if not query_seq or not target_seq:
            logging.warning(f"Empty sequence for {query_id}")
            results.append([query_id, target_contig, left_pos, right_pos, 3])
            continue
            
        target_insert_size = len(target_seq) - 200
        query_insert_size = len(query_seq) - 200
        if target_insert_size < query_insert_size * 0.3:
            presence = 1
        else:
            presence = 0

        results.append([query_id, target_contig, left_pos, right_pos, presence])
    
    return results

def main(args):
    # Parse BAM files for best alignments
    left_alignments = parse_best_alignments(args.left_bam)
    right_alignments = parse_best_alignments(args.right_bam)
    
    # Find ortholog pairs
    ortholog_pairs = find_ortholog_pairs(left_alignments, right_alignments)
    
    # Check locus presence
    results = check_locus_presence(args.query_fa, args.target_fa, ortholog_pairs)
    
    # Write results
    pd.DataFrame(results, columns=['ortholog_id', 'contig', 'left_flank_start', 
                                 'right_flank_end', 'presence']).astype({
        'left_flank_start': 'Int64',
        'right_flank_end': 'Int64',
        'presence': 'Int64'
    }).to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--left_bam')
    parser.add_argument('--right_bam')
    parser.add_argument('--query_fa')
    parser.add_argument('--target_fa')
    parser.add_argument('--output')
    args = parser.parse_args()
    
    main(args)