#!/usr/bin/env python3
"""
Generate per-read profiles from alignment based on variable positions.

This script:
1. Identifies variable positions in alignment
2. Creates per-read profiles based on those positions
3. Outputs profiles and metadata for downstream clustering

=====Script Structure=====
11 total functions organized into sections:

Utility Functions:
  1. log()                     - Centralized logging

File I/O Functions:
  2. read_fasta()              - BioPython FASTA reading
  3. seqs_to_array()           - Dict to numpy conversion
  4. array_to_seqs()           - Numpy to dict conversion
  5. validate_alignment()      - Validation logic

Processing Functions:
  6. mark_gap_extensions()     - Vectorized gap compression
  7. calculate_auto_trim()     - Vectorized trim calculation
  8. identify_variable_positions() - Optimized position finding
  9. create_read_profiles()    - Vectorized profile generation

Output Functions:
  10. log_parameters()         - Parameter logging
  11. write_profiles_tsv()     - Batch TSV writing

Main:
  12. main()                   - Clean workflow orchestration
"""

import sys
import argparse
from pathlib import Path
from collections import Counter
import numpy as np
from Bio import SeqIO


# ==================== Utility Functions ====================

def log(msg):
    """Simple stderr logger with newline."""
    sys.stderr.write(f"{msg}\n")


# ==================== File I/O Functions ====================

def read_fasta(path):
    """Read FASTA alignment using BioPython and return dict of {header: sequence}.
    
    Uses BioPython's optimized C-based parser for faster reading.
    Returns sequences as uppercase strings.
    """
    seqs = {}
    with open(path) as f:
        for record in SeqIO.parse(f, "fasta"):
            seqs[record.id] = str(record.seq).upper()
    return seqs


def seqs_to_array(seqs):
    """Convert sequence dictionary to numpy array and headers list.
    
    Returns:
        tuple: (headers_list, seq_array) where seq_array is 2D numpy array
    """
    headers = list(seqs.keys())
    seq_array = np.array([list(seq) for seq in seqs.values()], dtype='U1')
    return headers, seq_array


def array_to_seqs(headers, seq_array):
    """Convert numpy array back to sequence dictionary.
    
    Args:
        headers: List of sequence headers
        seq_array: 2D numpy array of sequences
    
    Returns:
        Dictionary of {header: sequence}
    """
    return {header: ''.join(seq_array[i]) for i, header in enumerate(headers)}


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
    log(f"Validated alignment: {len(seqs)} sequences, each {aln_len} bp")
    return aln_len


def mark_gap_extensions(seq_array, gap_char='-', extension_char='.'):
    """Prepare sequences with gap runs by marking continuation gaps with '.'.
    
    Operates directly on numpy array for efficiency.
    Replaces all gaps after the first in each consecutive gap run with '.'.
    This reduces the weight of long gaps in distance calculations while still
    penalizing the gap event itself.
    
    Example: 'A---T' becomes 'A-..T'
    
    Args:
        seq_array: 2D numpy array of sequences (rows=seqs, cols=positions)
        gap_char: Character representing gaps (default: '-')
        extension_char: Character to mark continuation gaps (default: '.')
    
    Returns:
        Modified 2D numpy array (operates in-place but also returns)
    """
    # Vectorized operations across all sequences at once
    is_gap = seq_array == gap_char
    
    # Find where gaps start (gap with no gap before)
    gap_starts = is_gap.copy()
    gap_starts[:, 1:] = is_gap[:, 1:] & ~is_gap[:, :-1]
    
    # Mark continuation gaps (gaps that aren't starts)
    continuation_gaps = is_gap & ~gap_starts
    
    # Replace continuation gaps with extension character
    seq_array[continuation_gaps] = extension_char
    
    return seq_array


def calculate_auto_trim(seq_array):
    """
    Calculate auto-trim values based on longest leading/trailing gaps.
    
    Vectorized implementation using numpy for efficiency.
    Returns the position of the first non-gap character in the sequence with
    the longest leading gap, and the position of the last non-gap character
    in the sequence with the longest trailing gap.
    
    Args:
        seq_array: 2D numpy array of sequences (rows=seqs, cols=positions)
    
    Returns:
        tuple: (trim_start, trim_end) - number of positions to trim from start and end
    """
    if seq_array.size == 0:
        return 0, 0
    
    # Find gap positions
    is_gap = seq_array == '-'
    
    # Find first non-gap position in each sequence
    # argmax finds first True in ~is_gap (first non-gap)
    has_non_gap = np.any(~is_gap, axis=1)
    first_non_gap = np.where(has_non_gap, np.argmax(~is_gap, axis=1), seq_array.shape[1])
    max_leading_gaps = int(first_non_gap.max())
    
    # Find last non-gap position in each sequence (search from right)
    # Reverse columns, find first non-gap, convert back to original position
    last_non_gap_from_end = np.where(has_non_gap, np.argmax(~is_gap[:, ::-1], axis=1), seq_array.shape[1])
    max_trailing_gaps = int(last_non_gap_from_end.max())
    
    return max_leading_gaps, max_trailing_gaps

def compute_trim_bounds(seq_array, auto_trim=True, trim_bp=70, min_trim=50, max_trim=250):
    """Compute start and end bounds (exclusive end) for variable position analysis.

    Centralizes all trim logic so downstream functions receive explicit bounds.

    Args:
        seq_array: 2D numpy array of sequences
        auto_trim: Use gap-derived automatic trimming if True
        trim_bp: Manual trim value (used if auto_trim False)
        min_trim: Minimum allowed trim (for auto mode)
        max_trim: Maximum allowed trim (for auto mode)

    Returns:
        tuple: (start_pos, end_pos, trim_start, trim_end, mode_str)
    """
    n_seqs, aln_len = seq_array.shape if seq_array.size else (0, 0)
    if seq_array.size == 0:
        return 0, 0, 0, 0, "empty"
    if auto_trim:
        raw_start, raw_end = calculate_auto_trim(seq_array)
        trim_start = max(min_trim, min(raw_start, max_trim))
        trim_end = max(min_trim, min(raw_end, max_trim))
        mode = "auto"
        log(f"Auto-trim detected: {raw_start} bp (start), {raw_end} bp (end)")
        log(f"Applied trim (constrained to [{min_trim}, {max_trim}]): {trim_start} bp (start), {trim_end} bp (end)")
    else:
        trim_start = trim_bp
        trim_end = trim_bp
        mode = "manual"
        log(f"Manual trim: ignoring first {trim_start} bp and last {trim_end} bp")
    start_pos = trim_start
    end_pos = max(start_pos, aln_len - trim_end)  # ensure non-negative region
    log(f"Analysis bounds: positions {start_pos}..{end_pos} (alignment length {aln_len})")
    return start_pos, end_pos, trim_start, trim_end, mode


def identify_variable_positions(seq_array, start_pos, end_pos, min_minor_freq, compress_gaps=True):
    """Identify positions whose second-most-common allele frequency >= threshold.

    Args:
        seq_array: 2D numpy array of sequences
        start_pos: Starting column index (inclusive)
        end_pos: Ending column index (exclusive)
        min_minor_freq: Threshold for second-most-common allele frequency
        compress_gaps: Skip positions where '.' is among top 2 alleles

    Returns:
        List of variable position indices (original alignment coordinates)
    """
    if seq_array.size == 0 or end_pos <= start_pos:
        return []
    n_seqs, _ = seq_array.shape
    region = seq_array[:, start_pos:end_pos]
    variable_positions = []
    # Iterate columns (vectorizing unique-with-count across all columns is non-trivial)
    for offset, pos in enumerate(range(start_pos, end_pos)):
        column = region[:, offset]
        unique, counts = np.unique(column, return_counts=True)
        if len(unique) < 2:
            continue
        if len(counts) == 2:
            sorted_idx = np.argsort(counts)[::-1]
            second_count = counts[sorted_idx[1]]
            top_2_chars = unique[sorted_idx[:2]]
        else:
            top_2_idx = np.argpartition(counts, -2)[-2:]
            sorted_top_2 = top_2_idx[np.argsort(counts[top_2_idx])[::-1]]
            second_count = counts[sorted_top_2[1]]
            top_2_chars = unique[sorted_top_2]
        if compress_gaps and '.' in top_2_chars:
            continue
        if (second_count / n_seqs) >= min_minor_freq:
            variable_positions.append(pos)
    return variable_positions


def create_read_profiles(seq_array, headers, variable_positions):
    """Create profile for each read based on variable positions.
    
    Works directly with numpy array for efficiency.
    
    Args:
        seq_array: 2D numpy array of sequences
        headers: List of sequence headers
        variable_positions: List of positions to extract
    
    Returns:
        Dictionary of {header: profile_tuple}
    """
    if not variable_positions:
        return {header: tuple() for header in headers}
    
    # Extract variable positions for all sequences at once (fancy indexing)
    var_pos_array = np.array(variable_positions)
    profiles_array = seq_array[:, var_pos_array]
    
    # Convert to dictionary of tuples
    profiles = {header: tuple(profiles_array[i]) for i, header in enumerate(headers)}
    
    return profiles


def log_parameters(args):
    """Log parameters to stderr."""
    log("=" * 60)
    log("Generate Read Profiles - Parameters")
    log("=" * 60)
    log(f"Input alignment: {args.alignment}")
    log(f"Output file: {args.output}")
    log(f"Min minor allele frequency: {args.min_minor_freq}")
    log(f"Trim mode: {'auto' if args.auto_trim else 'manual'}")
    if args.auto_trim:
        log(f"  Min trim: {args.min_trim} bp")
        log(f"  Max trim: {args.max_trim} bp")
    else:
        log(f"  Trim: {args.trim_bp} bp")
    log(f"Compress gaps: {args.compress_gaps}")
    log("=" * 60)
    log("")


def write_profiles_tsv(output_path, profiles, variable_positions):
    """Write profiles to TSV file with batch writing for performance.
    
    Args:
        output_path: Path object for output file
        profiles: Dictionary of {read_id: profile_tuple}
        variable_positions: List of variable position indices
    """
    log("Writing profiles to TSV...")
    
    # Build all lines in memory first
    lines = ["read_id\t" + "\t".join(f"pos_{p}" for p in variable_positions)]
    
    for read_id in sorted(profiles.keys()):
        profile = profiles[read_id]
        lines.append(f"{read_id}\t" + "\t".join(profile))
    
    # Single write operation
    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    
    log(f"Wrote {len(profiles)} read profiles to {output_path}")


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
    log_parameters(args)
    
    # Read and validate alignment
    log("Reading alignment...")
    seqs = read_fasta(args.alignment)
    log(f"Loaded {len(seqs)} sequences from alignment")
    
    try:
        aln_len = validate_alignment(seqs)
    except ValueError as e:
        log(f"ERROR: {e}")
        sys.exit(1)
    
    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Convert to numpy array once for efficiency
    log("\nConverting sequences to numpy array...")
    headers, seq_array = seqs_to_array(seqs)
    
    # Compress gap runs if requested
    if args.compress_gaps:
        log("Compressing gap runs (marking continuation gaps with '.')")
        seq_array = mark_gap_extensions(seq_array)
    
    # Compute trim bounds then identify variable positions
    log("\nComputing trim bounds...")
    start_pos, end_pos, trim_start, trim_end, trim_mode = compute_trim_bounds(
        seq_array,
        auto_trim=args.auto_trim,
        trim_bp=args.trim_bp,
        min_trim=args.min_trim,
        max_trim=args.max_trim,
    )
    log("\nIdentifying variable positions...")
    variable_positions = identify_variable_positions(
        seq_array, start_pos, end_pos, args.min_minor_freq, args.compress_gaps
    )
    log(f"Found {len(variable_positions)} variable positions (mode={trim_mode}, bounds {start_pos}..{end_pos}, threshold second allele >= {args.min_minor_freq})")
    
    # Warn about edge cases
    if len(variable_positions) == 0:
        log("WARNING: No variable positions found - all profiles will be identical")
    elif len(variable_positions) < 5:
        log(f"NOTE: Only {len(variable_positions)} variable positions found - "
            f"clustering may have limited resolution")
    
    # Create profiles
    log("\nGenerating per-read profiles...")
    profiles = create_read_profiles(seq_array, headers, variable_positions)
    
    # Write output
    log("")
    write_profiles_tsv(output_path, profiles, variable_positions)
    
    # Completion message
    log("\n" + "=" * 60)
    log("Profile generation complete")
    log("=" * 60)


if __name__ == "__main__":
    main()

