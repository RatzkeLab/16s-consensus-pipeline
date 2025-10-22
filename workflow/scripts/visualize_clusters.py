#!/usr/bin/env python3
"""
Generate per-cluster alignment files for visual inspection.

For each cluster identified by multi_consensus, this script extracts
the reads belonging to that cluster and saves them as a separate
alignment file.
"""

import sys
import argparse
from pathlib import Path
import re


def read_fasta(path):
    """Read FASTA file and return dict of {header: sequence}."""
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
                    seqs[header] = "".join(seq_buf)
                header = line[1:]
                seq_buf = []
            else:
                seq_buf.append(line)
        if header is not None:
            seqs[header] = "".join(seq_buf)
    return seqs


def parse_cluster_file(cluster_file):
    """
    Parse a cluster FASTA file to extract read IDs.
    
    The multi_consensus script names clusters like sample_A.fasta, sample_B.fasta
    We need to check if the file has a cluster_reads.txt or extract from the FASTA itself.
    """
    # For now, we'll determine cluster membership by looking for
    # corresponding _variants.tsv files which indicate multi-cluster output
    return None


def find_cluster_assignments(cluster_dir, sample_name):
    """
    Determine cluster assignments from multi_consensus output.
    
    Returns dict of {cluster_label: set_of_read_headers} or None if single cluster.
    """
    cluster_dir_path = Path(cluster_dir)
    
    # Find all cluster FASTA files
    cluster_files = list(cluster_dir_path.glob(f"{sample_name}_*.fasta"))
    
    if not cluster_files:
        # Single cluster - check for single file
        single_file = cluster_dir_path / f"{sample_name}.fasta"
        if single_file.exists():
            sys.stderr.write(f"Single cluster detected - no subclustering to visualize\n")
            return None
        else:
            sys.stderr.write(f"Warning: No cluster files found in {cluster_dir}\n")
            return None
    
    sys.stderr.write(f"Found {len(cluster_files)} cluster files\n")
    
    # For multi-cluster output, we need to infer which reads belong to which cluster
    # This requires re-running the clustering logic or reading a manifest file
    # For now, return cluster labels and we'll handle assignment separately
    clusters = {}
    for cluster_file in cluster_files:
        # Extract cluster label (e.g., "A" from "sample_A.fasta")
        match = re.match(rf"{sample_name}_([A-Z])", cluster_file.stem)
        if match:
            label = match.group(1)
            clusters[label] = cluster_file
    
    return clusters


def extract_cluster_assignments_from_log(log_file, sample_name):
    """
    Parse the multi_consensus log to extract cluster assignments.
    
    Log format includes lines like:
    "Cluster A: 47 reads -> A11_A.fasta"
    But doesn't include which specific reads - we need a different approach.
    """
    # This won't work without the actual read assignments
    # We need to modify multi_consensus.py to output a manifest file
    return None


def main():
    parser = argparse.ArgumentParser(
        description="Generate per-cluster alignments for visual inspection"
    )
    parser.add_argument("alignment", help="Input alignment FASTA file")
    parser.add_argument("cluster_dir", help="Multi-consensus output directory")
    parser.add_argument("outdir", help="Output directory for cluster alignments")
    parser.add_argument("sample", help="Sample name")
    
    args = parser.parse_args()
    
    # Read full alignment
    alignment = read_fasta(args.alignment)
    sys.stderr.write(f"Loaded {len(alignment)} sequences from alignment\n")
    
    # Check for cluster assignment file
    cluster_dir = Path(args.cluster_dir)
    assignment_file = cluster_dir / "cluster_assignments.tsv"
    
    if not assignment_file.exists():
        sys.stderr.write(f"No cluster assignments file found at {assignment_file}\n")
        sys.stderr.write(f"Multi-consensus likely produced single cluster - skipping visualization\n")
        
        # Create empty output directory
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        return
    
    # Read cluster assignments
    sys.stderr.write(f"Reading cluster assignments from {assignment_file}\n")
    cluster_assignments = {}  # read_header -> cluster_label
    
    with open(assignment_file) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                read_header = parts[0]
                cluster_label = parts[1]
                cluster_assignments[read_header] = cluster_label
    
    sys.stderr.write(f"Found assignments for {len(cluster_assignments)} reads\n")
    
    # Group reads by cluster
    clusters = {}
    for read_header, cluster_label in cluster_assignments.items():
        if cluster_label not in clusters:
            clusters[cluster_label] = []
        if read_header in alignment:
            clusters[cluster_label].append((read_header, alignment[read_header]))
    
    sys.stderr.write(f"Grouped reads into {len(clusters)} clusters\n")
    
    # Write per-cluster alignment files
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    for cluster_label, reads in sorted(clusters.items()):
        output_file = outdir / f"{args.sample}_cluster_{cluster_label}.fasta"
        
        with open(output_file, 'w') as f:
            for header, seq in reads:
                f.write(f">{header}\n")
                # Write sequence with line wrapping
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + "\n")
        
        sys.stderr.write(f"Wrote {len(reads)} reads to {output_file.name}\n")


if __name__ == "__main__":
    main()
