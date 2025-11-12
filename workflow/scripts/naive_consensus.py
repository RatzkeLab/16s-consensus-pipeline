#!/usr/bin/env python3
"""
Naive consensus sequence generation from multiple sequence alignment.

Generates a consensus sequence using a winner-takes-all approach (most common
base at each position) without ambiguous/wobble bases.

Usage:
    naive_consensus.py <alignment.fa> <output.fasta> <output.tsv> <sample_name> <min_prop>

Arguments:
    alignment.fa   - Input FASTA file with aligned sequences (equal length)
    output.fasta   - Output FASTA file for consensus sequence
    output.tsv     - Output TSV report of ambiguous positions
    sample_name    - Sample identifier for consensus sequence header
    min_prop       - Minimum proportion threshold (0-1) for calling a position
                     unambiguous. Positions where winner proportion <= this 
                     value are reported as ambiguous in the TSV.

Example:
    naive_consensus.py aligned.fa consensus.fa report.tsv sample_A 0.5
"""

import sys
import numpy as np

BASE_ORDER = ['A', 'T', 'G', 'C', '-']

def read_fasta(path):
    """Load FASTA sequences from file."""
    seqs = []
    with open(path) as f:
        seq_buf = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_buf:
                    seqs.append("".join(seq_buf).upper())
                seq_buf = []
            else:
                seq_buf.append(line)
        if seq_buf:
            seqs.append("".join(seq_buf).upper())
    return seqs

def build_frequency_matrix(seq_array):
    """Build base frequency matrix from sequence array."""
    n_seqs, aln_len = seq_array.shape
    freq_matrix = np.zeros((aln_len, len(BASE_ORDER)))
    
    for i, base in enumerate(BASE_ORDER):
        freq_matrix[:, i] = (seq_array == base).sum(axis=0) / n_seqs
    
    return freq_matrix

def get_consensus(freq_matrix):
    """Get consensus sequence and winner proportions from frequency matrix."""
    winner_indices = np.argmax(freq_matrix, axis=1)
    consensus = [BASE_ORDER[idx] for idx in winner_indices]
    winner_props = freq_matrix[np.arange(len(freq_matrix)), winner_indices]
    return consensus, winner_props

def find_ambiguous_positions(freq_matrix, winner_props, min_prop):
    """Find positions where winner proportion is below threshold."""
    ambiguous_mask = winner_props <= min_prop
    ambiguous_positions = []
    
    for pos in np.where(ambiguous_mask)[0]:
        variants = {base: freq_matrix[pos, i] 
                   for i, base in enumerate(BASE_ORDER) 
                   if freq_matrix[pos, i] > 0}
        ambiguous_positions.append((pos + 1, variants))
    
    return ambiguous_positions

def write_consensus_fasta(consensus, sample, output_path):
    """Write consensus sequence to FASTA file."""
    consensus_seq = "".join(consensus).replace("-", "")
    with open(output_path, "w") as f:
        f.write(f">{sample}\n")
        for i in range(0, len(consensus_seq), 80):
            f.write(consensus_seq[i:i+80] + "\n")

def write_ambiguous_report(ambiguous_positions, output_path):
    """Write ambiguous positions report to TSV file."""
    with open(output_path, "w") as f:
        f.write("position\tA\tT\tG\tC\tgap\n")
        for pos, variants in ambiguous_positions:
            a_prop = variants.get("A", 0.0)
            t_prop = variants.get("T", 0.0)
            g_prop = variants.get("G", 0.0)
            c_prop = variants.get("C", 0.0)
            gap_prop = variants.get("-", 0.0)
            f.write(f"{pos}\t{a_prop:.4f}\t{t_prop:.4f}\t{g_prop:.4f}\t{c_prop:.4f}\t{gap_prop:.4f}\n")

def validate_sequences(seqs, aln_path):
    """Validate that sequences are present and equal length."""
    if not seqs:
        sys.stderr.write(f"ERROR: No sequences in {aln_path}\n")
        sys.exit(1)
    
    aln_len = len(seqs[0])
    if any(len(s) != aln_len for s in seqs):
        sys.stderr.write("ERROR: Alignment sequences not equal length\n")
        sys.exit(1)
    
    return aln_len

def main(aln_path, fa_out, tsv_out, sample, min_prop):
    """Generate consensus sequence from alignment.
    
    Args:
        aln_path: Path to input FASTA file with aligned sequences
        fa_out: Path to output FASTA file for consensus sequence
        tsv_out: Path to output TSV file for ambiguous position report
        sample: Sample name for consensus sequence header
        min_prop: Minimum winner proportion (0-1) to call position unambiguous
    """
    seqs = read_fasta(aln_path)
    validate_sequences(seqs, aln_path)
    
    # Convert to numpy array
    seq_array = np.array([list(seq) for seq in seqs])
    
    # Build frequency matrix once
    freq_matrix = build_frequency_matrix(seq_array)
    
    # Get consensus and winner proportions
    consensus, winner_props = get_consensus(freq_matrix)
    
    # Find ambiguous positions
    ambiguous_positions = find_ambiguous_positions(freq_matrix, winner_props, min_prop)
    
    # Write outputs
    write_consensus_fasta(consensus, sample, fa_out)
    write_ambiguous_report(ambiguous_positions, tsv_out)

if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.stderr.write("Usage: naive_consensus.py <alignment> <out.fasta> <out.tsv> <sample> <min_prop>\n")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]))
