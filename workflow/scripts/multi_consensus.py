#!/usr/bin/env python3
"""
Multi-consensus generation for samples with multiple 16S copy numbers or contamination.

This script:
1. Identifies positions in alignment with low consensus (<min_agreement)
2. Creates per-read profile based on variants at these positions
3. Clusters reads using Hamming distance
4. Generates consensus for each cluster if min_cluster_size is met
"""

import sys
import argparse
from pathlib import Path
from collections import Counter, defaultdict
from itertools import combinations


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


def identify_variable_positions(seqs, min_agreement):
    """Identify positions where consensus is below threshold."""
    if not seqs:
        return []
    
    aln_len = len(list(seqs.values())[0])
    variable_positions = []
    
    for pos in range(aln_len):
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


def hamming_distance(profile1, profile2):
    """Calculate Hamming distance between two profiles."""
    return sum(c1 != c2 for c1, c2 in zip(profile1, profile2))


def cluster_reads(profiles, max_hamming):
    """Cluster reads based on Hamming distance using simple agglomerative approach."""
    # Start with each read in its own cluster
    clusters = {header: {header} for header in profiles.keys()}
    
    # Calculate all pairwise distances
    headers = list(profiles.keys())
    
    # Merge clusters if distance <= threshold
    merged = True
    while merged:
        merged = False
        cluster_ids = list(clusters.keys())
        
        for i, id1 in enumerate(cluster_ids):
            if id1 not in clusters:
                continue
            for id2 in cluster_ids[i+1:]:
                if id2 not in clusters:
                    continue
                
                # Check if any members of the two clusters are within threshold
                min_dist = float('inf')
                for h1 in clusters[id1]:
                    for h2 in clusters[id2]:
                        dist = hamming_distance(profiles[h1], profiles[h2])
                        min_dist = min(min_dist, dist)
                
                if min_dist <= max_hamming:
                    # Merge clusters
                    clusters[id1].update(clusters[id2])
                    del clusters[id2]
                    merged = True
                    break
            
            if merged:
                break
    
    return list(clusters.values())


def generate_consensus(seqs_dict, min_prop):
    """Generate consensus sequence from a subset of sequences."""
    seqs = list(seqs_dict.values())
    
    if not seqs:
        return "", []
    
    aln_len = len(seqs[0])
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
        
        consensus.append(winner_base)
        
        if winner_prop <= min_prop:
            variants = {base: count/total for base, count in counter.items()}
            ambiguous_positions.append((pos + 1, variants))
    
    # Remove gaps
    consensus_seq = "".join(consensus).replace("-", "")
    return consensus_seq, ambiguous_positions


def main():
    parser = argparse.ArgumentParser(description="Generate multi-consensus sequences")
    parser.add_argument("alignment", help="Input alignment FASTA file")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("sample", help="Sample name")
    parser.add_argument("--min_agreement", type=float, default=0.8)
    parser.add_argument("--min_cluster_size", type=int, default=5)
    parser.add_argument("--max_hamming", type=int, default=5)
    parser.add_argument("--min_consensus_prop", type=float, default=0.6)
    
    args = parser.parse_args()
    
    # Read alignment
    seqs = read_fasta(args.alignment)
    sys.stderr.write(f"Loaded {len(seqs)} sequences from alignment\n")
    
    if len(seqs) < args.min_cluster_size * 2:
        sys.stderr.write(f"Too few sequences for clustering. Creating single consensus.\n")
        consensus_seq, ambiguous = generate_consensus(seqs, args.min_consensus_prop)
        
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        
        with open(outdir / f"{args.sample}.fasta", "w") as f:
            f.write(f">{args.sample}\n")
            for i in range(0, len(consensus_seq), 80):
                f.write(consensus_seq[i:i+80] + "\n")
        
        sys.stderr.write(f"Single consensus written to {args.sample}.fasta\n")
        return
    
    # Identify variable positions
    variable_positions = identify_variable_positions(seqs, args.min_agreement)
    sys.stderr.write(f"Found {len(variable_positions)} variable positions\n")
    
    if len(variable_positions) < 3:
        sys.stderr.write(f"Too few variable positions. Creating single consensus.\n")
        consensus_seq, ambiguous = generate_consensus(seqs, args.min_consensus_prop)
        
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        
        with open(outdir / f"{args.sample}.fasta", "w") as f:
            f.write(f">{args.sample}\n")
            for i in range(0, len(consensus_seq), 80):
                f.write(consensus_seq[i:i+80] + "\n")
        
        sys.stderr.write(f"Single consensus written to {args.sample}.fasta\n")
        return
    
    # Create profiles
    profiles = create_read_profiles(seqs, variable_positions)
    
    # Cluster reads
    clusters = cluster_reads(profiles, args.max_hamming)
    sys.stderr.write(f"Found {len(clusters)} clusters\n")
    
    # Filter clusters by size
    valid_clusters = [c for c in clusters if len(c) >= args.min_cluster_size]
    sys.stderr.write(f"{len(valid_clusters)} clusters meet minimum size of {args.min_cluster_size}\n")
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if len(valid_clusters) < 2:
        sys.stderr.write(f"Only one valid cluster. Creating single consensus.\n")
        consensus_seq, ambiguous = generate_consensus(seqs, args.min_consensus_prop)
        
        with open(outdir / f"{args.sample}.fasta", "w") as f:
            f.write(f">{args.sample}\n")
            for i in range(0, len(consensus_seq), 80):
                f.write(consensus_seq[i:i+80] + "\n")
        
        sys.stderr.write(f"Single consensus written to {args.sample}.fasta\n")
        return
    
    # Generate consensus for each cluster
    cluster_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i, cluster in enumerate(sorted(valid_clusters, key=len, reverse=True)):
        label = cluster_labels[i] if i < len(cluster_labels) else str(i)
        cluster_seqs = {h: seqs[h] for h in cluster}
        
        consensus_seq, ambiguous = generate_consensus(cluster_seqs, args.min_consensus_prop)
        
        output_name = f"{args.sample}_{label}"
        with open(outdir / f"{output_name}.fasta", "w") as f:
            f.write(f">{output_name}\n")
            for j in range(0, len(consensus_seq), 80):
                f.write(consensus_seq[j:j+80] + "\n")
        
        # Write variants
        with open(outdir / f"{output_name}_variants.tsv", "w") as f:
            f.write("position\tA\tT\tG\tC\tgap\n")
            for pos, variants in ambiguous:
                a_prop = variants.get("A", 0.0)
                t_prop = variants.get("T", 0.0)
                g_prop = variants.get("G", 0.0)
                c_prop = variants.get("C", 0.0)
                gap_prop = variants.get("-", 0.0)
                f.write(f"{pos}\t{a_prop:.4f}\t{t_prop:.4f}\t{g_prop:.4f}\t{c_prop:.4f}\t{gap_prop:.4f}\n")
        
        sys.stderr.write(f"Cluster {label}: {len(cluster)} reads -> {output_name}.fasta\n")


if __name__ == "__main__":
    main()
