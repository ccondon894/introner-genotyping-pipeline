#!/usr/bin/env python3
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, Set, List, Tuple
import pandas as pd

samples = ["CCMP1545", "CCMP490", "RCC114", "RCC1614", "RCC1698", "RCC1749", "RCC2482", "RCC3052", "RCC373", "RCC465", "RCC629", "RCC647", "RCC692", "RCC693", "RCC833", "RCC835"]


def load_flank_pairs(left_fasta: Path, right_fasta: Path) -> Dict[str, str]:
    """Load corresponding left and right flank pairs from FASTA files.

    Args:
        left_fasta: Path to a FASTA file for left flanks
        right_fasta: Path to a FASTA file for right flanks

    Returns:
        flank_pairs: dict mapping left_header -> right_header
    """
    left_headers = []
    with open(left_fasta) as lf:
        for line in lf:
            if line.startswith('>'):
                left_headers.append(line.strip()[1:])  # remove '>'
                
    right_headers = []
    with open(right_fasta) as rf:
        for line in rf:
            if line.startswith('>'):
                right_headers.append(line.strip()[1:])
    
    # Validate equal number of sequences
    if len(left_headers) != len(right_headers):
        raise ValueError(
            f"Number of sequences in left fasta ({len(left_headers)}) "
            f"does not match right fasta ({len(right_headers)})"
        )
    
    # Create mapping of left -> right
    flank_pairs = dict(zip(left_headers, right_headers))
    return flank_pairs


def parse_seq_info(seq_id: str) -> Tuple[str, str, str, int, int, str]:
    split_pattern = r"#0#scaffold|-intronerized#0#|-deplete#0#|#0#intronerlesscontig|#0|#\d+"
    # Expected format: contig_seqid_L/R_start_end
    stuff = re.split(split_pattern, seq_id)
    sample_id = stuff[0]
    contig, the_rest = seq_id.strip().split('_introner_seq_')

    parts = the_rest.strip().split('_')
    seq_id = 'introner_seq_' + parts[0]
    flank_type = parts[1]
    start = int(parts[2])
    end = int(parts[3])
    if start > end:
        start, end = end, start

    return (sample_id, contig, flank_type, start, end, seq_id)




