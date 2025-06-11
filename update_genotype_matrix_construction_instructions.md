### Instructions for Adding Orientation Metadata to Genotype Matrix
Overview
We need to add orientation information to track whether left and right flanking sequences align to the forward strand or reverse complement during BWA alignment. This will help determine if flanks need to be flipped during pi/dxy calculations.
Implementation Plan
# 1. Modify the introner_genotype_matrix_builder.py script
A. Update the extract_best_alignments function

Add orientation detection logic to capture strand information from BAM alignments
Store orientation data alongside existing alignment information

Location: Around line 180-250 in the extract_best_alignments function
Changes needed:
python# In the alignment storage section, add orientation tracking
alignments[query_id] = {
    'query_coords': coords,
    'target_contig': bam.get_reference_name(read.reference_id),
    'target_start': read.reference_start,
    'target_end': read.reference_end,
    'quality': current_quality,
    'query_metadata': query_metadata,
    'is_reverse': read.is_reverse  # NEW: Add strand orientation
}
B. Update the build_ortholog_dictionary function

Modify the section where ortholog relationships are built (around line 300-400)
Capture orientation for both left and right flanking alignments
Determine overall orientation pattern for the introner

Logic to implement:
python# After extracting left and right alignments, determine orientation
left_reverse = left_aln.get('is_reverse', False) if left_aln else None
right_reverse = right_aln.get('is_reverse', False) if right_aln else None

# Determine overall orientation pattern
orientation = determine_orientation_pattern(left_reverse, right_reverse)
C. Add helper function for orientation determination
New function to add:
pythondef determine_orientation_pattern(left_reverse, right_reverse):
    """
    Determine the overall orientation pattern for an introner based on 
    left and right flank alignment orientations.
    
    Args:
        left_reverse: Boolean indicating if left flank aligned to reverse strand
        right_reverse: Boolean indicating if right flank aligned to reverse strand
    
    Returns:
        String indicating orientation pattern:
        - 'forward': Both flanks align forward
        - 'reverse': Both flanks align reverse (introner on reverse strand)
        - 'mixed': Mixed orientations (potential inversion/rearrangement)
        - 'partial': Only one flank available
        - 'unknown': Insufficient data
    """
    if left_reverse is None and right_reverse is None:
        return 'unknown'
    elif left_reverse is None or right_reverse is None:
        return 'partial'
    elif left_reverse == right_reverse:
        return 'reverse' if left_reverse else 'forward'
    else:
        return 'mixed'
D. Update the ortholog dictionary structure

Modify the dictionary that stores ortholog information (around line 450-500)
Add orientation fields to the data structure

Add to the ortholog_dict entries:
pythonortholog_dict[query][target][query_id] = {
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
    'target_metadata': target_metadata,
    'orientation': orientation,  # NEW: Add orientation information
    'left_reverse': left_reverse,  # NEW: Individual flank orientations
    'right_reverse': right_reverse  # NEW: Individual flank orientations
}
E. Update the output writing section

Modify the main output loop (around line 600-650)
Add new columns to the TSV header and data rows

Update header:
pythonoutfile.write("query\ttarget\tquery_id\tscenario\tquery_contig\tquery_start\tquery_end\ttarget_contig\t"
              "target_start\ttarget_end\tquery_gene\ttarget_gene\tquery_family\ttarget_family\t"
              "query_splice_site\ttarget_splice_site\tquery_introner_id\ttarget_introner_id\tmatch_score\t"
              "orientation\tleft_reverse\tright_reverse\n")  # NEW: Add orientation columns
Update data writing:
python# Extract orientation information
orientation = match.get('orientation', 'unknown')
left_reverse = match.get('left_reverse', '')
right_reverse = match.get('right_reverse', '')

# Convert boolean values to strings for output
left_reverse_str = str(left_reverse) if left_reverse is not None else ''
right_reverse_str = str(right_reverse) if right_reverse is not None else ''

# Add to the output line
outfile.write(f"{query}\t{target}\t{query_id}\t{match['scenario']}\t{match['query_contig'] or ''}\t"
              f"{query_start_str}\t{query_end_str}\t{match['target_contig'] or ''}\t"
              f"{target_start_str}\t{target_end_str}\t{query_gene}\t{target_gene}\t"
              f"{query_family}\t{target_family}\t{query_splice_site}\t{target_splice_site}\t"
              f"{query_introner_id}\t{target_introner_id}\t{match_score}\t"
              f"{orientation}\t{left_reverse_str}\t{right_reverse_str}\n")  # NEW: Add orientation data
# 2. Update the introner_network.py script
A. Modify the input processing

Update the script to handle the new orientation columns
Ensure orientation information is preserved through the network analysis

Location: Around line 150-200 in the processing loop
Changes needed:
python# Add orientation fields to the row processing
rows.append({
    'ortholog_id': ortholog_id,
    'sample': attrs['sample'],
    'start': attrs['start'],
    'end': attrs['end'],
    'contig': attrs['contig'],
    'sequence_id': node.split('|')[1] if '|' in node else '',
    'family': attrs.get('family', ''),
    'gene': attrs.get('gene', ''),
    'splice_site': attrs.get('splice_site', ''),
    'presence': attrs['presence'],
    'orientation': attrs.get('orientation', 'unknown'),  # NEW: Add orientation
    'left_reverse': attrs.get('left_reverse', ''),  # NEW: Add left flank orientation
    'right_reverse': attrs.get('right_reverse', '')  # NEW: Add right flank orientation
})

# Expected Outcome
After implementation, the genotype matrix will contain three new columns:

orientation: Overall orientation pattern ('forward', 'reverse', 'mixed', 'partial', 'unknown')
left_reverse: Boolean indicating if left flank aligned to reverse strand
right_reverse: Boolean indicating if right flank aligned to reverse strand

This information can then be used in pi/dxy calculations to determine when flanking sequences need to be reverse-complemented for proper analysis.