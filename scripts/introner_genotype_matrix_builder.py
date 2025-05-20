import os
import pysam
import pandas as pd
from collections import defaultdict
import sys

def load_introner_metadata(bed_dir, samples):
    """
    Load introner metadata from BED files.
    
    Returns: Dictionary mapping sample -> contig -> (start, end) -> metadata
    """
    metadata = defaultdict(lambda: defaultdict(dict))
    
    for sample in samples:
        bed_file = os.path.join(bed_dir, f"{sample}.candidate_loci_plus_flanks.filtered.similarity_checked.bed")
        if not os.path.exists(bed_file):
            print(f"Warning: BED file not found for {sample}")
            continue
            
        with open(bed_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 4:
                    continue
                    
                contig = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                introner_id = fields[3]
                
                # Parse metadata fields
                metadata_dict = {'introner_id': introner_id}
                for field in fields[4:]:
                    if '=' in field:
                        key, value = field.split('=', 1)
                        metadata_dict[key] = value
                
                # Store metadata using coordinates as key
                metadata[sample][contig][(start, end)] = metadata_dict
    
    return metadata

def find_matching_introner(metadata, sample, contig, position, wiggle=30):
    """
    Find an introner in the metadata that matches the given coordinates.
    
    Args:
        metadata: Introner metadata dictionary
        sample: Sample name
        contig: Contig name
        position: Position to check (start or end of alignment)
        wiggle: Wiggle room for coordinate matching in bp
        
    Returns:
        Tuple of (start, end) coordinates if found, None otherwise
    """
    if sample not in metadata or contig not in metadata[sample]:
        return None
        
    # Check each introner in this contig
    for (start, end), _ in metadata[sample][contig].items():
        # Check if position is near the left flank
        if abs(position - start) <= wiggle or abs(position - (start + 100)) <= wiggle:
            return (start, end)
            
        # Check if position is near the right flank
        if abs(position - (end - 100)) <= wiggle or abs(position - end) <= wiggle:
            return (start, end)
            
    return None

def find_overlapping_introner(metadata, sample, contig, target_start, target_end, min_overlap_fraction=0.3):
    """
    Find an introner in the metadata that overlaps with the target coordinates.
    
    Args:
        metadata: Introner metadata dictionary
        sample: Sample name
        contig: Contig name
        target_start: Start position of target region
        target_end: End position of target region
        min_overlap_fraction: Minimum fraction of overlap required (relative to the shorter region)
        
    Returns:
        Tuple of (start, end) coordinates of matching introner if found, None otherwise
    """
    if sample not in metadata or contig not in metadata[sample]:
        return None
    
    best_match = None
    best_overlap_fraction = 0
    
    # Check each introner in this contig
    for (start, end), _ in metadata[sample][contig].items():
        # Calculate overlap
        overlap_start = max(start, target_start)
        overlap_end = min(end, target_end)
        
        if overlap_start <= overlap_end:  # Regions overlap
            overlap_length = overlap_end - overlap_start + 1
            
            # Calculate overlap as a fraction of the shorter region
            target_length = target_end - target_start + 1
            introner_length = end - start + 1
            shorter_length = min(target_length, introner_length)
            
            overlap_fraction = overlap_length / shorter_length
            
            # Update best match if this one has a better overlap
            if overlap_fraction > best_overlap_fraction:
                best_overlap_fraction = overlap_fraction
                best_match = (start, end)
    
    # Return the best match if it meets the minimum overlap threshold
    if best_overlap_fraction >= min_overlap_fraction:
        return best_match
        
    # As a fallback, try using the original wiggle-based approach
    # This might catch cases where flanking alignments don't fully overlap with the introner
    wiggle = 35  # Same wiggle value as before
    for (start, end), _ in metadata[sample][contig].items():
        # Check if target_start is near the beginning of the introner
        if abs(target_start - start) <= wiggle:
            return (start, end)
        
        # Check if target_end is near the end of the introner
        if abs(target_end - end) <= wiggle:
            return (start, end)
            
    return None

def parse_query_name(query_name, flank_type):
    """Extract information from the query name."""
    # First, split by the pipe character
    parts = query_name.split('|')
    if len(parts) < 2:
        return None, None, {}
    
    # Get coordinates part (before the pipe)
    coords = parts[0]
    
    # Parse the coordinates to extract contig, start, end
    try:
        contig, pos_range = coords.split(':')
        start, end = map(int, pos_range.split('-'))
    except (ValueError, IndexError):
        contig, start, end = None, None, None
    
    # Get the introner ID and metadata (after the pipe)
    introner_info = parts[1].replace(f"_{flank_type}", "")
    
    # Split the introner info into ID and metadata
    introner_parts = introner_info.split()
    introner_id = introner_parts[0]
    
    # Parse metadata if available
    metadata = {}
    for part in introner_parts[1:]:
        if '=' in part:
            key, value = part.split('=', 1)
            metadata[key] = value
    
    # Create a unique ID by combining introner_id with coordinates
    unique_id = f"{introner_id}_{contig}_{start}_{end}" if contig and start and end else introner_id
    
    return coords, unique_id, metadata

def extract_best_alignments(bam_path, flank_type, export_quality_info=False):
    """
    Extract best alignments for each query from a BAM file.
    
    Args:
        bam_path: Path to the BAM file
        flank_type: Type of flank ('left' or 'right')
        export_quality_info: Whether to export detailed alignment quality info
        
    Returns: Dictionary mapping query_id -> alignment information
    """
    alignments = {}
    all_query_ids = set()  # Track all query IDs seen
    filtered_query_ids = set()  # Track query IDs retained after filtering
    low_quality_alignments = 0  # Track how many alignments fail quality checks
    alignment_quality_info = []  # Detailed quality info for problematic alignments
    
    if not os.path.exists(bam_path):
        return alignments
        
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
                
            coords, query_id, query_metadata = parse_query_name(read.query_name, flank_type)
            if not query_id:
                continue
            
            all_query_ids.add(query_id)
            
            # Check if alignment meets quality criteria
            if not is_good_alignment(read):
                low_quality_alignments += 1
                # If requested, collect detailed quality information for debugging
                if export_quality_info:
                    quality_info = analyze_alignment_quality(read)
                    quality_info['query_id'] = query_id
                    quality_info['flank_type'] = flank_type
                    alignment_quality_info.append(quality_info)
                continue
            
            # Store or update alignment if it's better than previous one
            current_quality = read.mapping_quality
            if query_id not in alignments or current_quality > alignments[query_id]['quality']:
                alignments[query_id] = {
                    'query_coords': coords,
                    'target_contig': bam.get_reference_name(read.reference_id),
                    'target_start': read.reference_start,
                    'target_end': read.reference_end,
                    'quality': current_quality,
                    'query_metadata': query_metadata
                }
                filtered_query_ids.add(query_id)
    
    # Print debug information
    print(f"BAM file: {bam_path}")
    print(f"Total unique query IDs found: {len(all_query_ids)}")
    print(f"Query IDs retained after quality filtering: {len(filtered_query_ids)}")
    print(f"Query IDs lost during filtering: {len(all_query_ids - filtered_query_ids)}")
    print(f"Low quality alignments rejected: {low_quality_alignments} (MAPQ < 30 or matched bases < 70)")
    
    if len(all_query_ids) > 0 and len(filtered_query_ids) > 0:
        # Print a sample of problematic query IDs (up to 5)
        lost_ids = list(all_query_ids - filtered_query_ids)
        if lost_ids:
            print(f"Sample of lost query IDs: {lost_ids[:min(5, len(lost_ids))]}")
    
    # Return different results depending on whether quality info was requested
    if export_quality_info:
        return alignments, alignment_quality_info
    else:
        return alignments

def build_ortholog_dictionary(samples, bam_dir, bed_dir, export_quality_info=False):
    """
    Build a comprehensive dictionary of orthologous relationships.
    
    Args:
        samples: List of sample names
        bam_dir: Directory containing BAM files
        bed_dir: Directory containing BED files
        export_quality_info: Whether to export detailed alignment quality info
    
    Returns: Dictionary mapping query -> target -> query_id -> alignment info
             And optionally a DataFrame with alignment quality information
    """
    # Load introner metadata from BED files
    print("Loading introner metadata...")
    introner_metadata = load_introner_metadata(bed_dir, samples)
    
    # Initialize result dictionary
    ortholog_dict = defaultdict(lambda: defaultdict(dict))
    
    # Track query IDs for debugging
    total_query_ids_processed = 0
    total_query_ids_kept = 0
    query_ids_by_scenario = {1: 0, 2: 0, 3: 0}  # Track counts by scenario
    
    # Track suspicious target regions
    suspicious_regions_count = 0
    
    # Collect alignment quality information if requested
    all_quality_info = []
    
    # Process each query-target pair
    print("Processing alignments...")
    for query in samples:
        for target in samples:
            if query == target:
                continue  # Skip self-comparisons
                
            print(f"Processing {query} -> {target}")
            
            # Get left and right flank alignments
            left_bam = os.path.join(bam_dir, f"{query}_{target}.left_flank.sorted.bam")
            right_bam = os.path.join(bam_dir, f"{query}_{target}.right_flank.sorted.bam")
            if not os.path.exists(left_bam) or not os.path.exists(right_bam):
                print(f"Warning: BAM files not found for {query} -> {target}")
                print(left_bam)
                print(right_bam)
                continue  # Skip this pair if files don't exist
                
            # Extract alignments, potentially with quality info
            if export_quality_info:
                left_alignments, left_quality_info = extract_best_alignments(left_bam, "left", export_quality_info=True)
                right_alignments, right_quality_info = extract_best_alignments(right_bam, "right", export_quality_info=True)
                
                # Collect quality info with query/target identifiers
                for info in left_quality_info + right_quality_info:
                    info['query_sample'] = query
                    info['target_sample'] = target
                    all_quality_info.append(info)
            else:
                left_alignments = extract_best_alignments(left_bam, "left")
                right_alignments = extract_best_alignments(right_bam, "right")
            
            # Find query IDs in either left or right alignments
            all_query_ids = set(left_alignments.keys()) | set(right_alignments.keys())
            print(f"Total unique query IDs for {query} -> {target}: {len(all_query_ids)}")
            
            # Count query IDs with only left or right alignment
            left_only = set(left_alignments.keys()) - set(right_alignments.keys())
            right_only = set(right_alignments.keys()) - set(left_alignments.keys())
            both = set(left_alignments.keys()) & set(right_alignments.keys())
            print(f"Left flank only: {len(left_only)}, Right flank only: {len(right_only)}, Both flanks: {len(both)}")
            
            total_query_ids_processed += len(all_query_ids)
            
            # Process each query ID
            for query_id in all_query_ids:
                left_aln = left_alignments.get(query_id)
                right_aln = right_alignments.get(query_id)
                
                # Check for the specific problematic query mentioned
                is_problematic_query = "CCMP1545#0#scaffold_15:648456-648819|introner_seq_4078" in query_id
                
                # Get query metadata
                query_contig, coords_range = left_aln['query_coords'].split(':') if left_aln else (right_aln['query_coords'].split(':') if right_aln else (None, None))
                query_start, query_end = map(int, coords_range.split('-')) if coords_range else (None, None)
                query_coords = None
                
                if not left_aln or not right_aln:
                    ortholog_dict[query][target][query_id] = {
                        'scenario': 3,
                        'query_id': query_id,
                        'query_coords': left_aln['query_coords'] if left_aln else (right_aln['query_coords'] if right_aln else None),
                        'query_contig': query_contig,
                        'query_start': query_start,
                        'query_end': query_end,
                        'target_contig': None,
                        'target_start': None,
                        'target_end': None,
                        'query_metadata': None,
                        'target_metadata': None
                    }
                    query_ids_by_scenario[3] += 1
                    continue
                
                # Check if alignments are on same contig but not in correct orientation
                if left_aln['target_contig'] != right_aln['target_contig']:
                    ortholog_dict[query][target][query_id] = {
                        'scenario': 3,
                        'query_id': query_id,
                        'query_coords': left_aln['query_coords'],
                        'query_contig': query_contig,
                        'query_start': query_start,
                        'query_end': query_end,
                        'target_contig': None,
                        'target_start': None,
                        'target_end': None,
                        'query_metadata': None,
                        'target_metadata': None
                    }
                    query_ids_by_scenario[3] += 1
                    continue
                    
                # If alignments are on the same contig but inverted (left end > right start),
                # we'll still process them as they could be due to genome rearrangements
                is_inverted = left_aln['target_end'] > right_aln['target_start']
                
                # Case 1 or 2: Both flanks aligned (potentially inverted)
                target_contig = left_aln['target_contig']
                
                # For inverted alignments, flip the coordinates for consistency
                if is_inverted:
                    target_start = right_aln['target_start']
                    target_end = left_aln['target_end']
                else:
                    target_start = left_aln['target_start']
                    target_end = right_aln['target_end']
                
                # Debug output for the problematic query or any suspicious ranges
                target_length = (target_end - target_start) if target_end and target_start else 0
                query_length = query_end - query_start if query_end and query_start else 0
                
                # Define size thresholds for flagging suspicious pairs
                min_target_length = 50  # Minimum acceptable target length
                min_ratio = 0.2  # Target should be at least this fraction of query length
                
                if (is_problematic_query or 
                    (target_length < min_target_length and query_length > 5 * target_length)):
                    suspicious_regions_count += 1
                    print(f"\nSUSPICIOUS REGION #{suspicious_regions_count}: {query} -> {target}, {query_id}")
                    print(f"  Query: {query_contig}:{query_start}-{query_end} (length: {query_length})")
                    print(f"  Left alignment: {left_aln['target_contig']}:{left_aln['target_start']}-{left_aln['target_end']}")
                    print(f"  Right alignment: {right_aln['target_contig']}:{right_aln['target_start']}-{right_aln['target_end']}")
                    print(f"  Calculated target: {target_contig}:{target_start}-{target_end} (length: {target_length})")
                    print(f"  Is inverted: {is_inverted}")
                    
                    # For very small targets, try an alternative approach - expand to include the full region between alignments
                    if target_length < min_target_length:
                        if is_inverted:
                            # For inverted alignments, take the full span
                            corrected_start = min(right_aln['target_start'], left_aln['target_start'])
                            corrected_end = max(right_aln['target_end'], left_aln['target_end'])
                        else:
                            # For normal orientation, take region between the alignments
                            corrected_start = left_aln['target_end']
                            corrected_end = right_aln['target_start']
                        
                        # Apply correction if it actually helps
                        corrected_length = corrected_end - corrected_start
                        
                        # Only use correction if it makes the target larger and more proportional to query
                        if corrected_length > target_length and corrected_length >= min(query_length * min_ratio, min_target_length):
                            print(f"  CORRECTION APPLIED: New target: {target_contig}:{corrected_start}-{corrected_end} (new length: {corrected_length})")
                            target_start = corrected_start
                            target_end = corrected_end
                        else:
                            # If correction doesn't help, try to ensure minimum reasonable target size
                            # based on query size, if we have target metadata
                            mid_point = (target_start + target_end) // 2
                            desired_length = max(min_target_length, int(query_length * min_ratio))
                            half_length = desired_length // 2
                            
                            # Only apply this expansion if we have confidence in the region
                            if target_length > 0:
                                target_start = mid_point - half_length
                                target_end = mid_point + half_length
                                print(f"  SIZE CORRECTION: New target: {target_contig}:{target_start}-{target_end} (new length: {target_end-target_start})")
                
                for (start, end), meta in introner_metadata[query][query_contig].items():
                    if start <= query_start <= end and start <= query_end <= end:
                        query_coords = (start, end)
                        break
                
                # Get query metadata - prefer from BAM query name, fall back to BED file
                query_name_metadata = left_aln.get('query_metadata', {}) or right_aln.get('query_metadata', {})
                
                # Check if target has an introner at this location using the overlap-based approach
                target_coords = find_overlapping_introner(
                    introner_metadata, target, target_contig, target_start, target_end
                )
                
                # If no overlapping introner found, fall back to the position-based approach
                if target_coords is None:
                    target_coords = find_matching_introner(
                        introner_metadata, target, target_contig, target_start
                    )
                
                # Get metadata from BED files
                query_bed_metadata = introner_metadata[query][query_contig].get(query_coords) if query_coords else None
                target_metadata = introner_metadata[target][target_contig].get(target_coords) if target_coords else None
                
                # Combine the metadata, prioritizing the query name metadata
                query_metadata = {}
                if query_bed_metadata:
                    query_metadata.update(query_bed_metadata)
                if query_name_metadata:
                    query_metadata.update(query_name_metadata)
                
                # Determine scenario
                scenario = 1 if target_metadata else 2
                
                # Use target coordinates from metadata if found, otherwise use alignment coordinates
                if target_coords and scenario == 1:
                    target_start_final = target_coords[0]
                    target_end_final = target_coords[1]
                    
                    # Debug output for the case where metadata coords are used
                    if is_problematic_query:
                        print(f"  Using metadata coordinates: {target_contig}:{target_start_final}-{target_end_final} (length: {target_end_final-target_start_final})")
                else:
                    target_start_final = target_start
                    target_end_final = target_end
                
                ortholog_dict[query][target][query_id] = {
                    'scenario': scenario,
                    'query_id': query_id,
                    'query_coords': left_aln['query_coords'],
                    'query_contig': query_contig,
                    'query_start': query_start,
                    'query_end': query_end,
                    'target_contig': target_contig,
                    'target_start': target_start_final,
                    'target_end': target_end_final,
                    'query_metadata': query_metadata,
                    'target_metadata': target_metadata
                }
                query_ids_by_scenario[scenario] += 1
                total_query_ids_kept += 1
            
            # Add unaligned introners from metadata
            ortholog_dict = add_unaligned_introners_to_ortholog_dict(ortholog_dict, introner_metadata, query, target)
    
    # Print debug information
    print(f"Total query IDs processed: {total_query_ids_processed}")
    print(f"Total query IDs kept: {total_query_ids_kept}")
    print(f"Suspicious regions detected and corrected: {suspicious_regions_count}")
    for scenario, count in query_ids_by_scenario.items():
        print(f"Scenario {scenario}: {count} query IDs")
    
    # Return different results depending on whether quality info was requested
    if export_quality_info and all_quality_info:
        quality_df = pd.DataFrame(all_quality_info)
        return ortholog_dict, quality_df
    else:
        return ortholog_dict

def select_best_matches(ortholog_dict):
    """
    For cases with multiple potential targets, select the best matches based on metadata.
    """
    for query in ortholog_dict:
        for target in ortholog_dict[query]:
            for query_id, match in ortholog_dict[query][target].items():
                # Skip cases without both metadata
                if match['scenario'] != 1 or not match['query_metadata'] or not match['target_metadata']:
                    continue
                    
                # Compare gene and introner_id
                score = 0
                
                # Check gene match
                q_gene = match['query_metadata'].get('gene')
                t_gene = match['target_metadata'].get('gene')
                if q_gene and t_gene and q_gene == t_gene:
                    score += 10
                    
                # Check introner_id match
                q_id = match['query_metadata'].get('introner_id')
                t_id = match['target_metadata'].get('introner_id')
                if q_id and t_id and q_id == t_id:
                    score += 5
                    
                # Store score
                match['match_score'] = score
                
    return ortholog_dict

def get_all_introners_from_metadata(metadata, query):
    """
    Extract all introners for a given query sample from metadata.
    
    Args:
        metadata: Dictionary of introner metadata
        query: Query sample name
        
    Returns:
        List of dictionaries with introner information
    """
    introners = []
    
    if query not in metadata:
        return introners
        
    # Iterate through all contigs and introners for this query
    for contig, introner_dict in metadata[query].items():
        for (start, end), meta_dict in introner_dict.items():
            # Create a unique ID for this introner
            introner_id = meta_dict.get('introner_id', 'unknown')
            unique_id = f"{introner_id}_{contig}_{start}_{end}"
            
            # Create introner entry
            introner_info = {
                'query_id': unique_id,
                'query_coords': f"{contig}:{start}-{end}",
                'query_contig': contig,
                'query_start': start,
                'query_end': end,
                'query_metadata': meta_dict
            }
            
            introners.append(introner_info)
            
    return introners

def add_unaligned_introners_to_ortholog_dict(ortholog_dict, introner_metadata, query, target):
    """
    Add introners from metadata that don't have alignments to the ortholog dictionary.
    
    Args:
        ortholog_dict: Dictionary of aligned introners
        introner_metadata: Dictionary of introner metadata
        query: Query sample name
        target: Target sample name
        
    Returns:
        Updated ortholog dictionary
    """
    # Get all introners for this query from metadata
    all_introners = get_all_introners_from_metadata(introner_metadata, query)
    
    # Get set of query_ids already in the ortholog dictionary for this pair
    existing_query_ids = set(ortholog_dict[query][target].keys())
    
    # Count of added introners
    added_count = 0
    
    # Add introners not yet in ortholog_dict
    for introner in all_introners:
        query_id = introner['query_id']
        
        if query_id not in existing_query_ids:
            # Add as scenario 3 (no alignment)
            ortholog_dict[query][target][query_id] = {
                'scenario': 3,
                'query_id': query_id,
                'query_coords': introner['query_coords'],
                'query_contig': introner['query_contig'],
                'query_start': introner['query_start'],
                'query_end': introner['query_end'],
                'target_contig': None,
                'target_start': None,
                'target_end': None,
                'query_metadata': introner['query_metadata'],
                'target_metadata': None
            }
            added_count += 1
            
    if added_count > 0:
        print(f"  Added {added_count} unaligned introners from metadata for {query} -> {target}")
        
    return ortholog_dict

def is_good_alignment(read, min_mapq=30, min_match_length=70):
    """
    Check if an alignment is of sufficient quality.
    
    Args:
        read: pysam AlignedSegment object
        min_mapq: Minimum mapping quality score (default: 30)
        min_match_length: Minimum number of matched bases in CIGAR (default: 70)
    
    Returns:
        Boolean indicating if alignment meets quality criteria
    """
    # Check mapping quality
    if read.mapping_quality < min_mapq:
        return False
    
    # Calculate total matched bases from CIGAR
    matched_bases = 0
    for op, length in read.cigartuples:
        # CIGAR operations: 0=M (match/mismatch), 7== (match), 8=X (mismatch)
        if op in [0, 7, 8]:
            matched_bases += length
    
    # Check if matched bases meet minimum threshold
    return matched_bases >= min_match_length

def analyze_alignment_quality(read):
    """
    Analyze an alignment's quality and return diagnostic information.
    
    Args:
        read: pysam AlignedSegment object
    
    Returns:
        Dictionary with alignment quality metrics
    """
    # CIGAR operation types
    cigar_ops = {
        0: 'M',  # match or mismatch
        1: 'I',  # insertion
        2: 'D',  # deletion
        3: 'N',  # skipped region
        4: 'S',  # soft clip
        5: 'H',  # hard clip
        6: 'P',  # padding
        7: '=',  # sequence match
        8: 'X'   # sequence mismatch
    }
    
    # Count each operation type
    op_counts = defaultdict(int)
    total_length = 0
    
    for op, length in read.cigartuples:
        op_char = cigar_ops.get(op, '?')
        op_counts[op_char] += length
        total_length += length if op != 5 else 0  # exclude hard clips from total length
    
    # Calculate aligned length (M/=/X operations)
    aligned_length = op_counts.get('M', 0) + op_counts.get('=', 0) + op_counts.get('X', 0)
    
    # Calculate alignment percentage (aligned / total sequence length)
    align_percent = (aligned_length / total_length) * 100 if total_length > 0 else 0
    
    return {
        'mapq': read.mapping_quality,
        'cigar_str': read.cigarstring,
        'aligned_length': aligned_length,
        'total_length': total_length,
        'align_percent': f"{align_percent:.1f}%",
        'op_counts': dict(op_counts),
        'read_name': read.query_name,
        'reverse_strand': read.is_reverse,
        'ref_name': read.reference_name,
        'ref_start': read.reference_start,
        'ref_end': read.reference_end
    }

def main():
    # Configuration
    samples = ["CCMP1545", "RCC114", "RCC1614", "RCC1698", "RCC1749", 
               "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", 
               "RCC692", "RCC693", "RCC833"]
    
    bam_dir = sys.argv[1]
    bed_dir = sys.argv[1]
    
    # Check if we should export alignment quality info
    export_quality_info = False
    quality_output_file = None
    if len(sys.argv) > 3 and sys.argv[3] == "--export-quality":
        export_quality_info = True
        quality_output_file = sys.argv[4] if len(sys.argv) > 4 else "alignment_quality_info.tsv"
        print(f"Will export alignment quality info to {quality_output_file}")
    
    # Build the ortholog dictionary
    if export_quality_info:
        ortholog_dict, quality_df = build_ortholog_dictionary(samples, bam_dir, bed_dir, export_quality_info=True)
    else:
        ortholog_dict = build_ortholog_dictionary(samples, bam_dir, bed_dir)
    
    # Select best matches based on metadata
    ortholog_dict = select_best_matches(ortholog_dict)
    
    # Count entries before writing
    total_entries = 0
    query_entries_by_sample = {sample: 0 for sample in samples}
    
    for query in ortholog_dict:
        for target in ortholog_dict[query]:
            entries_count = len(ortholog_dict[query][target])
            total_entries += entries_count
            query_entries_by_sample[query] = query_entries_by_sample.get(query, 0) + entries_count
    
    print(f"Total entries to write: {total_entries}")
    print("Entries by query sample:")
    for sample, count in query_entries_by_sample.items():
        print(f"  {sample}: {count}")
    
    # Output results to file
    with open(sys.argv[2], "w") as outfile:
        # Write header
        outfile.write("query\ttarget\tquery_id\tscenario\tquery_contig\tquery_start\tquery_end\ttarget_contig\t"
                      "target_start\ttarget_end\tquery_gene\ttarget_gene\tquery_family\ttarget_family\t"
                      "query_splice_site\ttarget_splice_site\tquery_introner_id\ttarget_introner_id\tmatch_score\n")
        
        # Track rows written
        rows_written = 0
        
        # Write data
        for query in ortholog_dict:
            for target in ortholog_dict[query]:
                for query_id, match in ortholog_dict[query][target].items():
                    # Extract metadata fields with empty defaults if not available
                    query_gene = match['query_metadata'].get('gene', '') if match['query_metadata'] else ''
                    target_gene = match['target_metadata'].get('gene', '') if match['target_metadata'] else ''
                    
                    # Extract family values and ensure they're integers if present
                    query_family = match['query_metadata'].get('family', '') if match['query_metadata'] else ''
                    if query_family and str(query_family).strip() and query_family != '':
                        try:
                            query_family = int(float(query_family))
                        except (ValueError, TypeError):
                            pass
                            
                    target_family = match['target_metadata'].get('family', '') if match['target_metadata'] else ''
                    if target_family and str(target_family).strip() and target_family != '':
                        try:
                            target_family = int(float(target_family))
                        except (ValueError, TypeError):
                            pass
                    
                    query_splice_site = match['query_metadata'].get('splice_site', '') if match['query_metadata'] else ''
                    target_splice_site = match['target_metadata'].get('splice_site', '') if match['target_metadata'] else ''
                    query_introner_id = match['query_metadata'].get('introner_id', '') if match['query_metadata'] else ''
                    target_introner_id = match['target_metadata'].get('introner_id', '') if match['target_metadata'] else ''
                    match_score = match.get('match_score', 0)
                    
                    # Format start and end values as integers
                    query_start = match['query_start']
                    query_end = match['query_end']
                    target_start = match['target_start']
                    target_end = match['target_end']
                    
                    # Ensure start is always smaller than end for both query and target
                    if query_start is not None and query_end is not None and query_start > query_end:
                        query_start, query_end = query_end, query_start
                        
                    if target_start is not None and target_end is not None and target_start > target_end:
                        target_start, target_end = target_end, target_start
                    
                    # Convert numeric values to integer format for output
                    query_start_str = str(int(query_start)) if query_start is not None else ''
                    query_end_str = str(int(query_end)) if query_end is not None else ''
                    target_start_str = str(int(target_start)) if target_start is not None else ''
                    target_end_str = str(int(target_end)) if target_end is not None else ''
                    
                    outfile.write(f"{query}\t{target}\t{query_id}\t{match['scenario']}\t{match['query_contig'] or ''}\t"
                                  f"{query_start_str}\t{query_end_str}\t{match['target_contig'] or ''}\t"
                                  f"{target_start_str}\t{target_end_str}\t{query_gene}\t{target_gene}\t"
                                  f"{query_family}\t{target_family}\t{query_splice_site}\t{target_splice_site}\t"
                                  f"{query_introner_id}\t{target_introner_id}\t{match_score}\n")
                    rows_written += 1
    
    print(f"Total rows written to output file: {rows_written}")
    print("Done! Results written to ortholog_results.tsv")
    
    # Export alignment quality info if available
    if export_quality_info and 'quality_df' in locals():
        quality_df.to_csv(quality_output_file, sep='\t', index=False)
        print(f"Exported alignment quality information to {quality_output_file}")
    
    return ortholog_dict

if __name__ == "__main__":
    main()