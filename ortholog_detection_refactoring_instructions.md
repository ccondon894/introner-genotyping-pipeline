# Ortholog Detection Pipeline Refactoring Instructions

## Completed Tasks
- **Phase 1**: Implemented new parsing functions in `introner_genotype_matrix_builder.py`
  - Added `parse_whole_genome_alignment()` to analyze SAM files from whole genome alignments
  - Added `classify_introner_alignment()` to determine introner presence/absence based on alignment patterns
  - Added `build_ortholog_dictionary_whole_genome()` to process all samples using the new method
  - Updated `build_ortholog_dictionary()` to support whole-genome mode
  - Added command-line flag `--whole-genome-mode` using argparse
  - Enhanced output format with additional columns for alignment quality, coverage, and gap size

- **Phase 2**: Updated the Snakemake workflow in `ortholog_detection_pipeline.smk`
  - Added whole-genome alignment rules using minimap2
  - Created `pair_orthologs_from_whole_genome` rule to process all-vs-all alignments
  - Added `build_genotype_matrix_whole_genome` rule to generate the whole-genome based matrix
  - Defined sample groups (Group 1 and Group 2) for appropriate alignment parameters
  - Set up correct alignment presets (asm5/asm20) based on within/between group comparisons

## Pending Tasks
- **Phase 3**: Testing the new implementation
  - Test with a subset of 3-4 samples to validate the approach
  - Compare results from old and new methods 
  - Manually inspect a subset of classifications
  - Benchmark performance (runtime and memory usage)

- **Phase 4**: Final validation and cleanup
  - Confirm that scenario classifications are more consistent than current method
  - Ensure gene-based scoring provides meaningful ortholog groupings
  - Verify missing data rates are similar or lower than current approach
  - Remove old flank-based code once new method is confirmed working

## Overview
This document outlines the refactoring of `scripts/introner_genotype_matrix_builder.py` and `ortholog_detection_pipeline.smk` to use whole genome alignments instead of separate flank alignments for improved ortholog detection accuracy. The new approach will perform all-vs-all genome alignments and parse alignment results for putative introner regions to determine presence/absence patterns.

## Current Workflow Structure
- Current workflow extracts left and right flanks from candidate introner loci
- Performs separate BWA alignments of flanks to target genomes
- Infers ortholog relationships by analyzing flank alignment positions
- Builds ortholog networks through complex coordinate matching logic

## Proposed Changes

### 1. Replace Flank-Based Alignment with Whole Genome Alignment

#### Current Implementation Issues:
- Separate flank alignments can lead to coordinate mismatches
- Limited context around introner sites
- Complex logic for inferring presence/absence from flank positions
- Difficulty handling structural variants and repetitive regions

#### New Implementation:
- Perform all-vs-all whole genome alignments using minimap2
- Parse SAM files directly for introner regions of interest
- Use alignment gaps/presence to determine ortholog status

### 2. Minimap2 Alignment Strategy

#### Alignment Parameters:
```bash
# For internal group sample pairs (Group 2 vs Group 2 and Group 1 vs Group 1)
minimap2 -x asm5 -a --eqx -t {threads} {reference}.fa {query}.fa > {output}.sam

# For between Group 2 vs. Group 1 comparisons (RCC1749 and RCC3052 vs. Group 1)
minimap2 -x asm20 -a --eqx -t {threads} {reference}.fa {query}.fa > {output}.sam
```

#### Sample Group Definitions:
- **Group 1**: CCMP1545, RCC114, RCC1614, RCC1698, RCC2482, RCC373, RCC465, RCC629, RCC692, RCC693, RCC833
- **Group 2**: RCC1749, RCC3052

#### Alignment Logic:
- Use `asm5` when query and target are fro mthe same group
- Use `asm20` when both query and target are from different groups

### 3. SAM File Parsing for Introner Regions

#### Region Definition:
For each candidate introner locus, define the region of interest as:
- **Query region**: introner coordinates + 100bp flanks on each side
- **Target region**: corresponding aligned region in target genome

#### Alignment Classification Logic:
1. **Present in both (Scenario 1)**: 
   - Continuous alignment across entire region (flanks + introner)
   - No large gaps or insertions in target alignment
   
2. **Absent in target (Scenario 2)**:
   - Flanks align well to target genome
   - Gap/deletion in target genome between aligned flanks
   - Gap size approximately matches introner length in query
   
3. **Missing data (Scenario 3)**:
   - Poor or no alignment of flanking regions
   - Insufficient alignment quality to make determination
   - Target region not covered in alignment

### 4. Implementation Details for introner_genotype_matrix_builder.py

#### New Function Structure:
```python
def parse_whole_genome_alignment(sam_file, query_introner_coords, min_mapq=30):
    """
    Parse SAM file for introner region alignment patterns
    
    Args:
        sam_file: Path to minimap2 SAM output
        query_introner_coords: Dict of {introner_id: (contig, start, end)}
        min_mapq: Minimum mapping quality threshold
    
    Returns:
        Dict of alignment classifications for each introner
    """
    
def classify_introner_alignment(alignment_blocks, query_start, query_end):
    """
    Classify introner presence/absence based on alignment pattern
    
    Args:
        alignment_blocks: List of aligned segments from SAM
        query_start: Start coordinate of introner + 100bp flank
        query_end: End coordinate of introner + 100bp flank
    
    Returns:
        scenario: 1 (present), 2 (absent), 3 (missing)
        target_coords: Coordinates in target genome if applicable
    """
```

#### SAM Parsing Logic:
1. **Read SAM file and filter by MAPQ ≥ 30**
2. **For each introner region of interest:**
   - Extract all alignment blocks overlapping the region (introner + 200bp total)
   - Analyze alignment pattern:
     - **Continuous coverage**: Scenario 1 (present in both)
     - **Gap in middle with flanks covered**: Scenario 2 (absent in target)
     - **Poor/no coverage**: Scenario 3 (missing data)

#### Coordinate Handling:
- **Query coordinates**: Use introner coordinates from metadata + 100bp flanks
- **Target coordinates**: Extract from SAM alignment coordinates
- **Gap detection**: Look for deletions/unaligned regions between well-aligned flanks

### 5. Gene-Based Scoring Enhancement

#### Maintain Current Gene Matching Logic:
```python
# Enhanced scoring system (keep existing logic)
score = 0

# Check gene match (existing logic)
q_gene = match['query_metadata'].get('gene')
t_gene = match['target_metadata'].get('gene')
if q_gene and t_gene and q_gene == t_gene:
    score += 10
    
# Check introner_id match (existing logic)  
q_id = match['query_metadata'].get('introner_id')
t_id = match['target_metadata'].get('introner_id')
if q_id and t_id and q_id == t_id:
    score += 5

# Store score
match['match_score'] = score
```

### 6. Updated Workflow Integration

#### New Rule Structure for ortholog_detection_pipeline.smk:
```python
rule whole_genome_alignment:
    input:
        query_genome = os.path.join(GENOME_DIR, "{query}.vg_paths.fa"),
        target_genome = os.path.join(GENOME_DIR, "{target}.vg_paths.fa")
    output:
        sam = os.path.join(OUTPUT_DIR, "{query}_{target}.whole_genome.sam")
    params:
        preset = lambda wildcards: "asm20" if (wildcards.query in ["RCC1749", "RCC3052"] and wildcards.target in ["RCC1749", "RCC3052"]) else "asm5"
    threads: 8
    shell:
        """
        minimap2 -x {params.preset} -a --eqx -t {threads} {input.target_genome} {input.query_genome} > {output.sam}
        """

rule pair_orthologs_from_whole_genome:
    input:
        sam_files = [os.path.join(OUTPUT_DIR, f"{query}_{target}.whole_genome.sam")
                     for query in SAMPLES for target in SAMPLES if query != target],
        candidate_beds = expand(os.path.join(OUTPUT_DIR, "{sample}.candidate_loci_plus_flanks.filtered.similarity_checked.bed"), sample=SAMPLES)
    output:
        orthologs = os.path.join(GT_DIR, "ortholog_results.tsv")
    shell:
        """
        python scripts/introner_genotype_matrix_builder.py {OUTPUT_DIR} {output.orthologs} --whole-genome-mode
        """
```

### 7. Error Handling and Edge Cases

#### Handle Common Issues:
- **Multiple alignments**: Choose best alignment by MAPQ score
- **Repetitive regions**: Use alignment uniqueness metrics
- **Structural variants**: Detect complex rearrangements that may affect classification
- **Assembly gaps**: Mark as missing data when target regions contain N's

#### Quality Control:
- Log alignment statistics for each sample pair
- Track alignment coverage across introner regions
- Report problematic regions that may need manual inspection

### 8. Output Format Compatibility

#### Maintain Existing Output Structure:
Keep the same TSV output format with columns:
- query, target, query_id, scenario, query_contig, query_start, query_end
- target_contig, target_start, target_end, query_gene, target_gene
- query_family, target_family, query_splice_site, target_splice_site
- query_introner_id, target_introner_id, match_score

#### Enhanced Logging:
Add additional columns for debugging:
- alignment_quality: MAPQ score of best alignment
- alignment_coverage: Fraction of introner region covered
- gap_size: Size of gap in target genome (for scenario 2)

## Implementation Notes

### Performance Considerations:
- Whole genome alignments will take longer than flank alignments
- Consider parallelizing alignment jobs across sample pairs
- SAM files will be larger but can be processed streaming

### Testing Strategy:
1. **Test with subset**: Start with 3-4 samples to validate approach
2. **Compare results**: Run both old and new methods on test data
3. **Validate scenarios**: Manually inspect a subset of classifications
4. **Benchmark performance**: Measure runtime and memory usage

### Validation Criteria:
- Scenario classifications should be more consistent than current method
- Gene-based scoring should still provide meaningful ortholog groupings
- Missing data rates should be similar or lower than current approach

## Migration Path

1. **Phase 1**: Implement new parsing functions in introner_genotype_matrix_builder.py
2. **Phase 2**: Add command-line flag for whole-genome mode
3. **Phase 3**: Update Snakemake rules to use new alignment strategy
4. **Phase 4**: Validate results and remove old flank-based code

This refactoring should provide more accurate ortholog detection while maintaining compatibility with downstream analysis steps.