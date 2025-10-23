#!/usr/bin/env python3
"""
Split subsampled reads into cluster-specific FASTQ files.

Reads cluster assignments and splits the subsampled FASTQ file
into separate files for each cluster.
"""

import sys
import argparse
from pathlib import Path


def read_cluster_assignments(assignment_file):
    """Read cluster assignments from TSV file."""
    assignments = {}
    with open(assignment_file) as f:
        header = f.readline()  # skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                # Extract just the UUID (first part before space) to match FASTQ parsing
                full_id = parts[0]
                read_id = full_id.split()[0] if ' ' in full_id else full_id
                cluster = parts[1]
                assignments[read_id] = cluster
    return assignments


def parse_fastq(fastq_path):
    """Generator that yields (header, sequence, quality) tuples from FASTQ."""
    with open(fastq_path) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            
            # Extract read ID (everything after @ and before first space)
            read_id = header[1:].split()[0] if header.startswith('@') else header
            
            yield read_id, header, seq, qual


def main():
    parser = argparse.ArgumentParser(description="Split reads by cluster assignment")
    parser.add_argument("fastq", help="Input subsampled FASTQ file")
    parser.add_argument("cluster_dir", help="Cluster detection output directory")
    parser.add_argument("outdir", help="Output directory for split FASTQs")
    parser.add_argument("sample", help="Sample name")
    
    args = parser.parse_args()
    
    cluster_dir = Path(args.cluster_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Check if clustering was performed
    no_clusters_file = cluster_dir / "no_clusters.txt"
    assignment_file = cluster_dir / "cluster_assignments.tsv"
    
    if no_clusters_file.exists():
        # No clustering - just copy the input file
        sys.stderr.write("No clusters detected - copying input file\n")
        import shutil
        shutil.copy(args.fastq, outdir / f"{args.sample}.fastq")
        sys.stderr.write(f"Copied to {args.sample}.fastq\n")
        return
    
    if not assignment_file.exists():
        sys.stderr.write("ERROR: No cluster assignments or no_clusters marker found\n")
        sys.exit(1)
    
    # Read cluster assignments
    assignments = read_cluster_assignments(assignment_file)
    sys.stderr.write(f"Loaded assignments for {len(assignments)} reads\n")
    
    # Open output file handles for each cluster
    cluster_files = {}
    clusters_seen = set(assignments.values())
    
    for cluster in sorted(clusters_seen):
        output_path = outdir / f"{args.sample}_{cluster}.fastq"
        cluster_files[cluster] = open(output_path, 'w')
        sys.stderr.write(f"Writing cluster {cluster} to {output_path.name}\n")
    
    # Split reads
    reads_written = {cluster: 0 for cluster in clusters_seen}
    reads_unassigned = 0
    
    for read_id, header, seq, qual in parse_fastq(args.fastq):
        if read_id in assignments:
            cluster = assignments[read_id]
            cluster_files[cluster].write(f"{header}\n{seq}\n+\n{qual}\n")
            reads_written[cluster] += 1
        else:
            reads_unassigned += 1
    
    # Close files
    for f in cluster_files.values():
        f.close()
    
    # Report
    sys.stderr.write(f"\nSplit complete:\n")
    for cluster in sorted(clusters_seen):
        sys.stderr.write(f"  Cluster {cluster}: {reads_written[cluster]} reads\n")
    
    if reads_unassigned > 0:
        sys.stderr.write(f"  Unassigned: {reads_unassigned} reads (not in alignment)\n")


if __name__ == "__main__":
    main()
