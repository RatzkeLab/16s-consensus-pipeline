#!/usr/bin/env python3
"""
Generate per-read profiles from alignment based on variable positions.

This script:
1. Identifies variable positions in alignment
2. Creates per-read profiles based on those positions
3. Outputs profiles and metadata for downstream clustering
"""

import sys
import argparse
from pathlib import Path
from collections import Counter
import numpy as np

# Import common helper
sys.path.insert(0, str(Path(__file__).parent))
from common_helpers import calculate_auto_trim


def read_fasta(path):
    """Read FASTA alignment and return dict of {header: sequence}.
    
    Optimized to build sequences character-by-character with list.extend()
    and single join at the end per sequence.
    """
    seqs = {}
    with open(path) as f:
        header = None
        seq_chars = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = ''.join(seq_chars)
                header = line[1:]
                seq_chars = []
            else:
                # Extend with uppercase characters (no intermediate strings)
                seq_chars.extend(line.upper())
        if header is not None:
            seqs[header] = ''.join(seq_chars)
    return seqs


def validate_alignment(seqs):
    """Validate that alignment is non-empty and all sequences have same length.
    
    Args:
        seqs: Dictionary of {header: sequence}
    
    Raises:
        ValueError: If alignment is empty or sequences have inconsistent lengths
    """
    if not seqs:
        raise ValueError("Alignment file is empty or contains no sequences")
    
    lengths = {len(seq) for seq in seqs.values()}
    
    if len(lengths) > 1:
        # Report which sequences have different lengths
        length_info = {}
        for header, seq in seqs.items():
            seq_len = len(seq)
            if seq_len not in length_info:
                length_info[seq_len] = []
            length_info[seq_len].append(header)
        
        error_msg = "Alignment sequences have inconsistent lengths:\n"
        for length, headers in sorted(length_info.items()):
            error_msg += f"  Length {length}: {len(headers)} sequences (e.g., {headers[0]})\n"
        raise ValueError(error_msg)
    
    aln_len = lengths.pop()
    sys.stderr.write(f"Validated alignment: {len(seqs)} sequences, each {aln_len} bp\n")
    return aln_len


def mark_gap_extensions(seqs, gap_char='-', extension_char='.'):
    """Prepare sequences with gap runs by marking continuation gaps with '.'.
    
    Replaces all gaps after the first in each consecutive gap run with '.'.
    This reduces the weight of long gaps in distance calculations while still
    penalizing the gap event itself.
    
    Example: 'A---T' becomes 'A-..T'
    
    Args:
        seqs: Dictionary of {header: sequence}
        gap_char: Character representing gaps (default: '-')
        extension_char: Character to mark continuation gaps (default: '.')
    
    Returns:
        Dictionary of {header: modified_sequence}
    """
    compressed = {}
    
    for header, seq in seqs.items():
        # Convert to numpy array for vectorized operations
        seq_array = np.array(list(seq), dtype='U1')
        
        # Find gap positions
        is_gap = seq_array == gap_char
        
        # Find where gaps start (gap with no gap before)
        gap_starts = is_gap.copy()
        gap_starts[1:] = is_gap[1:] & ~is_gap[:-1]
        
        # Mark continuation gaps (gaps that aren't starts)
        continuation_gaps = is_gap & ~gap_starts
        
        # Replace continuation gaps with extension character
        seq_array[continuation_gaps] = extension_char
        
        compressed[header] = ''.join(seq_array)
    
    return compressed


def identify_variable_positions(seqs, min_minor_freq, trim_bp=70, auto_trim=True, compress_gaps=True, 
                               min_trim=50, max_trim=250):
    """Identify positions where the second-most-common allele is above threshold.
    
    Uses numpy for vectorized operations on large alignments.
    
    Args:
        seqs: Dictionary of sequences
        min_minor_freq: Minimum frequency for second-most-common allele (default: 0.05)
        trim_bp: Number of bp to ignore at start and end of alignment (default: 70)
        auto_trim: If True, automatically calculate trim values based on alignment gaps
        compress_gaps: If True, filter out '.' (gap continuation marker) when it's in top 2 most common
        min_trim: Minimum trim value - auto_trim cannot trim less than this (default: 50)
        max_trim: Maximum trim value - auto_trim cannot trim more than this (default: 250)
    """
    if not seqs:
        return []
    
    # Convert sequences to numpy array for vectorized operations
    seq_list = list(seqs.values())
    aln_len = len(seq_list[0])
    n_seqs = len(seq_list)
    
    # Create 2D array: rows = sequences, cols = positions
    # Use Unicode strings with max length 1
    seq_array = np.array([list(seq) for seq in seq_list], dtype='U1')
    
    variable_positions = []
    
    # Calculate trim values
    if auto_trim:
        auto_trim_start, auto_trim_end = calculate_auto_trim(seqs)
        # Constrain auto-trim to be between min_trim and max_trim
        trim_start = max(min_trim, min(auto_trim_start, max_trim))
        trim_end = max(min_trim, min(auto_trim_end, max_trim))
        sys.stderr.write(f"Auto-trim detected: {auto_trim_start} bp (start), {auto_trim_end} bp (end)\n")
        sys.stderr.write(f"Applied trim (constrained to [{min_trim}, {max_trim}]): "
                        f"{trim_start} bp (start), {trim_end} bp (end)\n")
    else:
        trim_start = trim_bp
        trim_end = trim_bp
        sys.stderr.write(f"Manual trim: ignoring first {trim_start} bp and last {trim_end} bp\n")
    
    # Define the region to analyze (excluding trimmed regions)
    start_pos = trim_start
    end_pos = aln_len - trim_end
    
    sys.stderr.write(f"Analyzing positions {start_pos} to {end_pos}\n")

    # Iterate over each position in the alignment, and save variable positions
    for pos in range(start_pos, end_pos):
        column = seq_array[:, pos]
        
        # Get unique characters and their counts
        unique, counts = np.unique(column, return_counts=True)

        if len(unique) < 2: 
            # Only one character type at this position
            continue

        # Sort by count (descending)
        sorted_indices = np.argsort(counts)[::-1]
        sorted_counts = counts[sorted_indices]
        
        # Check second-most-common frequency
        second_count = sorted_counts[1]
        second_freq = second_count / len(column)

        # If compress_gaps is enabled, check if '.' is in top 2 most common
        # Filter it out if so
        if compress_gaps:
            if '.' in unique:
                # Get top 2 by count
                top_2_chars = unique[sorted_indices[:2]]
                if '.' in top_2_chars:
                    # skip this position, as gap continuation is prevalent
                    continue
        
        if second_freq >= min_minor_freq:
            variable_positions.append(pos)
    
    return variable_positions


def create_read_profiles(seqs, variable_positions):
    """Create profile for each read based on variable positions."""
    profiles = {}
    for header, seq in seqs.items():
        profile = tuple(seq[pos] for pos in variable_positions)
        profiles[header] = profile
    return profiles


def main():
    parser = argparse.ArgumentParser(description="Generate per-read profiles from alignment")
    parser.add_argument("alignment", help="Input alignment FASTA file")
    parser.add_argument("output", help="Output TSV file")
    parser.add_argument("--min_minor_freq", type=float, default=0.05,
                        help="Minimum frequency for second-most-common allele to mark position as variable (default: 0.05)")
    parser.add_argument("--trim_bp", type=int, default=70,
                        help="Number of bp to ignore at start and end of alignment")
    parser.add_argument("--auto_trim", action="store_true",
                        help="Automatically calculate trim values based on longest leading/trailing gaps")
    parser.add_argument("--min_trim", type=int, default=50,
                        help="Minimum trim value - auto_trim cannot trim less than this (default: 50)")
    parser.add_argument("--max_trim", type=int, default=250,
                        help="Maximum trim value - auto_trim cannot trim more than this (default: 250)")
    parser.add_argument("--compress_gaps", action="store_true",
                        help="Compress gap runs by marking continuation gaps (reduces weight of long gaps)")
    
    args = parser.parse_args()
    
    # Log parameters
    sys.stderr.write("=" * 60 + "\n")
    sys.stderr.write("Generate Read Profiles - Parameters\n")
    sys.stderr.write("=" * 60 + "\n")
    sys.stderr.write(f"Input alignment: {args.alignment}\n")
    sys.stderr.write(f"Output file: {args.output}\n")
    sys.stderr.write(f"Min minor allele frequency: {args.min_minor_freq}\n")
    sys.stderr.write(f"Trim mode: {'auto' if args.auto_trim else 'manual'}\n")
    if args.auto_trim:
        sys.stderr.write(f"  Min trim: {args.min_trim} bp\n")
        sys.stderr.write(f"  Max trim: {args.max_trim} bp\n")
    else:
        sys.stderr.write(f"  Trim: {args.trim_bp} bp\n")
    sys.stderr.write(f"Compress gaps: {args.compress_gaps}\n")
    sys.stderr.write("=" * 60 + "\n\n")
    
    # Read alignment
    sys.stderr.write("Reading alignment...\n")
    seqs = read_fasta(args.alignment)
    sys.stderr.write(f"Loaded {len(seqs)} sequences from alignment\n")
    
    # Validate alignment
    try:
        aln_len = validate_alignment(seqs)
    except ValueError as e:
        sys.stderr.write(f"ERROR: {e}\n")
        sys.exit(1)
    
    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Compress gap runs if requested
    if args.compress_gaps:
        sys.stderr.write("\nCompressing gap runs (marking continuation gaps with '.')\n")
        seqs = mark_gap_extensions(seqs)
    
    # Identify variable positions
    sys.stderr.write("\nIdentifying variable positions...\n")
    variable_positions = identify_variable_positions(
        seqs, args.min_minor_freq, args.trim_bp, args.auto_trim, args.compress_gaps,
        args.min_trim, args.max_trim
    )
    sys.stderr.write(f"Found {len(variable_positions)} variable positions "
                    f"(threshold: second allele >= {args.min_minor_freq})\n")
    
    if len(variable_positions) == 0:
        sys.stderr.write("WARNING: No variable positions found - all profiles will be identical\n")
    elif len(variable_positions) < 5:
        sys.stderr.write(f"NOTE: Only {len(variable_positions)} variable positions found - "
                        f"clustering may have limited resolution\n")
    
    # Create profiles (even if empty - let clustering script decide what to do)
    sys.stderr.write("\nGenerating per-read profiles...\n")
    profiles = create_read_profiles(seqs, variable_positions)
    
    # Write read profiles to TSV
    with open(output_path, "w") as f:
        # Header: read_id, then each variable position
        f.write("read_id\t" + "\t".join(f"pos_{p}" for p in variable_positions) + "\n")
        
        # Each read's profile
        for read_id in sorted(profiles.keys()):
            profile = profiles[read_id]
            f.write(f"{read_id}\t" + "\t".join(profile) + "\n")
    
    sys.stderr.write(f"Wrote {len(profiles)} read profiles to {output_path}\n")
    sys.stderr.write("\n" + "=" * 60 + "\n")
    sys.stderr.write("Profile generation complete\n")
    sys.stderr.write("=" * 60 + "\n")


if __name__ == "__main__":
    main()

