#!/usr/bin/env python3
"""
Calculate pairwise edit distances between consensus sequences.

Reads a multi-FASTA file and computes the edit distance (Levenshtein distance)
between every pair of sequences, outputting a TSV matrix.
"""

import sys
from pathlib import Path


def parse_fasta(fasta_path):
    """Parse FASTA file and return dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Add last sequence
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)
    
    return sequences


def edit_distance(s1, s2):
    """Calculate Levenshtein edit distance between two strings."""
    m, n = len(s1), len(s2)
    
    # Create DP table
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    # Initialize base cases
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j
    
    # Fill DP table
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s1[i-1] == s2[j-1]:
                dp[i][j] = dp[i-1][j-1]
            else:
                dp[i][j] = 1 + min(
                    dp[i-1][j],    # deletion
                    dp[i][j-1],    # insertion
                    dp[i-1][j-1]   # substitution
                )
    
    return dp[m][n]


def main():
    # Get input/output from snakemake
    input_fasta = snakemake.input.multi_db
    output_tsv = snakemake.output.distances
    
    # Parse sequences
    sequences = parse_fasta(input_fasta)
    seq_names = list(sequences.keys())
    
    # Calculate pairwise distances
    with open(output_tsv, 'w') as out:
        # Write header
        out.write("Sequence1\tSequence2\tEditDistance\n")
        
        # Calculate all pairwise distances
        for i, name1 in enumerate(seq_names):
            for j, name2 in enumerate(seq_names):
                if i < j:  # Only calculate upper triangle (symmetric matrix)
                    dist = edit_distance(sequences[name1], sequences[name2])
                    out.write(f"{name1}\t{name2}\t{dist}\n")
    
    print(f"Calculated {len(seq_names) * (len(seq_names) - 1) // 2} pairwise distances", 
          file=sys.stderr)


if __name__ == "__main__":
    main()
