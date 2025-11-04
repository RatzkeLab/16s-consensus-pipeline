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

# Import common helper
sys.path.insert(0, str(Path(__file__).parent))
from common_helpers import calculate_auto_trim


def read_fasta(path):
    """Read FASTA alignment and return dict of {header: sequence}."""
    seqs = {}
    with open(path) as f:
        header = None
        seq_buf = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs[header] = "".join(seq_buf).upper()
                header = line[1:]
                seq_buf = []
            else:
                seq_buf.append(line)
        if header is not None:
            seqs[header] = "".join(seq_buf).upper()
    return seqs


def compress_gap_runs(seqs):
    """Compress gap runs by marking continuation gaps with 'N'.
    
    Replaces all gaps after the first in each consecutive gap run with 'N'.
    This reduces the weight of long gaps in distance calculations while still
    penalizing the gap event itself.
    
    Example: 'A---T' becomes 'A-NNT'
    
    Args:
        seqs: Dictionary of {header: sequence}
    
    Returns:
        Dictionary of {header: modified_sequence}
    """
    compressed = {}
    
    for header, seq in seqs.items():
        new_seq = []
        in_gap_run = False
        
        for char in seq:
            if char == '-':
                if in_gap_run:
                    # This is a continuation gap, mark with N
                    new_seq.append('N')
                else:
                    # This is the first gap in a run, keep it
                    new_seq.append('-')
                    in_gap_run = True
            else:
                # Non-gap character, reset gap run flag
                new_seq.append(char)
                in_gap_run = False
        
        compressed[header] = ''.join(new_seq)
    
    return compressed


def identify_variable_positions(seqs, min_agreement, trim_bp=70, auto_trim=False):
    """Identify positions where consensus is below threshold.
    
    Args:
        seqs: Dictionary of sequences
        min_agreement: Maximum agreement threshold for variable positions
        trim_bp: Number of bp to ignore at start and end of alignment (default: 70)
        auto_trim: If True, automatically calculate trim values based on alignment gaps
    """
    if not seqs:
        return []
    
    aln_len = len(list(seqs.values())[0])
    variable_positions = []
    
    # Calculate trim values
    if auto_trim:
        trim_start, trim_end = calculate_auto_trim(seqs)
        sys.stderr.write(f"Auto-trim enabled: ignoring first {trim_start} bp and last {trim_end} bp "
                        f"(based on longest leading/trailing gaps)\n")
    else:
        trim_start = trim_bp
        trim_end = trim_bp
        sys.stderr.write(f"Manual trim: ignoring first {trim_start} bp and last {trim_end} bp\n")
    
    # Define the region to analyze (excluding trimmed regions)
    start_pos = trim_start
    end_pos = aln_len - trim_end
    
    sys.stderr.write(f"Analyzing positions {start_pos} to {end_pos}\n")
    
    for pos in range(start_pos, end_pos):
        column = [seq[pos] for seq in seqs.values()]
        counter = Counter(column)
        total = len(column)
        
        most_common_count = counter.most_common(1)[0][1]
        agreement = most_common_count / total
        
        if agreement <= min_agreement:
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
    parser.add_argument("--min_agreement", type=float, default=0.8,
                        help="Maximum agreement threshold for identifying variable positions")
    parser.add_argument("--trim_bp", type=int, default=70,
                        help="Number of bp to ignore at start and end of alignment")
    parser.add_argument("--auto_trim", action="store_true",
                        help="Automatically calculate trim values based on longest leading/trailing gaps")
    parser.add_argument("--compress_gaps", action="store_true",
                        help="Compress gap runs by marking continuation gaps (reduces weight of long gaps)")
    
    args = parser.parse_args()
    
    # Read alignment
    seqs = read_fasta(args.alignment)
    sys.stderr.write(f"Loaded {len(seqs)} sequences from alignment\n")
    
    # Ensure output directory exists
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Compress gap runs if requested
    if args.compress_gaps:
        sys.stderr.write("Compressing gap runs (marking continuation gaps with 'N')\n")
        seqs = compress_gap_runs(seqs)
    
    # Identify variable positions
    variable_positions = identify_variable_positions(
        seqs, args.min_agreement, args.trim_bp, args.auto_trim
    )
    sys.stderr.write(f"Found {len(variable_positions)} variable positions\n")
    
    # Create profiles (even if empty - let clustering script decide what to do)
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
    sys.stderr.write("Profile generation complete\n")


if __name__ == "__main__":
    main()
