#!/usr/bin/env python3
"""
Naive consensus sequence generation from multiple sequence alignment.
Winner-takes-all approach without wobble bases.
"""

import sys
from collections import Counter

def read_fasta(path):
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

def main(aln_path, fa_out, tsv_out, sample, min_prop):
    seqs = read_fasta(aln_path)
    
    if not seqs:
        sys.stderr.write(f"ERROR: No sequences in {aln_path}\n")
        sys.exit(1)
    
    aln_len = len(seqs[0])
    if any(len(s) != aln_len for s in seqs):
        sys.stderr.write("ERROR: Alignment sequences not equal length\n")
        sys.exit(1)
    
    consensus = []
    ambiguous_positions = []
    
    for pos in range(aln_len):
        column = [seq[pos] for seq in seqs]
        counter = Counter(column)
        total = len(column)
        
        most_common = counter.most_common(1)[0]
        winner_base = most_common[0]
        winner_count = most_common[1]
        winner_prop = winner_count / total
        
        # Always use winner base (no N's), but record ambiguous positions
        consensus.append(winner_base)
        
        if winner_prop <= min_prop:
            variants = {base: count/total for base, count in counter.items()}
            ambiguous_positions.append((pos + 1, variants))
    
    consensus_seq = "".join(consensus).replace("-", "")
    with open(fa_out, "w") as f:
        f.write(f">{sample}\n")
        for i in range(0, len(consensus_seq), 80):
            f.write(consensus_seq[i:i+80] + "\n")
    
    with open(tsv_out, "w") as f:
        f.write("position\tA\tT\tG\tC\tgap\n")
        for pos, variants in ambiguous_positions:
            a_prop = variants.get("A", 0.0)
            t_prop = variants.get("T", 0.0)
            g_prop = variants.get("G", 0.0)
            c_prop = variants.get("C", 0.0)
            gap_prop = variants.get("-", 0.0)
            f.write(f"{pos}\t{a_prop:.4f}\t{t_prop:.4f}\t{g_prop:.4f}\t{c_prop:.4f}\t{gap_prop:.4f}\n")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.stderr.write("Usage: naive_consensus.py <alignment> <out.fasta> <out.tsv> <sample> <min_prop>\n")
        sys.exit(1)
    
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]))
