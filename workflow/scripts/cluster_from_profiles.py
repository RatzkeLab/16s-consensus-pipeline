#!/usr/bin/env python3
"""
Detect if a sample has obvious subclusters based on pre-generated read profiles.

This script:
1. Reads per-read profiles from upstream step
2. Performs hierarchical clustering
3. Determines if multiple valid clusters exist
4. Outputs cluster assignments IF clusters are found
"""

"""Cluster detection from per-read profiles.

Imports are organized into standard library, third-party scientific stack.
"""

# Standard library
import sys
import argparse
from pathlib import Path
from collections import defaultdict

# Third-party
import numpy as np
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend for pipeline/server environments
from scipy.cluster.hierarchy import linkage, fcluster

# Local visualization helpers separated out
from cluster_visualization import (
    plot_distance_heatmap,
    plot_profiles_clustermap,
    write_trivial_profiles_viz,
)

# ==================== Utility ====================

def log(msg: str) -> None:
    """Unified stderr logging (single place to adjust formatting)."""
    sys.stderr.write(f"{msg}\n")


def read_profiles(profile_file):
    """Read read profiles from TSV file.
    
    Returns:
        Tuple of (profiles_dict, variable_positions)
        - profiles_dict: {read_id: tuple of characters at variable positions}
        - variable_positions: list of position numbers (parsed from header)
    """
    profiles = {}
    variable_positions = []
    
    with open(profile_file) as f:
        # Parse header to get variable positions
        header = f.readline().strip().split('\t')
        variable_positions = [int(pos.replace('pos_', '')) for pos in header[1:]]
        
        # Parse each read's profile
        for line in f:
            if not line.strip():
                continue
            fields = line.strip().split('\t')
            read_id = fields[0]
            profile = tuple(fields[1:])
            profiles[read_id] = profile
    
    return profiles, variable_positions


def hamming_distance(profile1, profile2):
    """Calculate Hamming-like distance between two profiles.

    Neutral characters ('.','N') are ignored; mismatches among remaining positions count as 1.
    """
    dist = 0
    for c1, c2 in zip(profile1, profile2):
        if c1 in {'.','N'} or c2 in {'.','N'}:
            continue
        if c1 != c2:
            dist += 1
    return dist

def compute_distance_matrix(profiles, neutral_chars={'.','N'}, max_tensor_cells=2_50_000_000):
    """Vectorized pairwise distance computation with memory guard."""
    headers = list(profiles.keys())
    n = len(headers)
    if n == 0:
        return np.zeros((0,0))
    m = len(profiles[headers[0]]) if profiles[headers[0]] else 0
    if n == 1 or m == 0:
        return np.zeros((n,n))
    if n * n * m > max_tensor_cells:
        dist_matrix = np.zeros((n,n), dtype=float)
        for i in range(n):
            pi = profiles[headers[i]]
            for j in range(i+1, n):
                pj = profiles[headers[j]]
                dist = hamming_distance(pi, pj)
                dist_matrix[i,j] = dist
                dist_matrix[j,i] = dist
        return dist_matrix
    char_array = np.array([profiles[h] for h in headers], dtype='U1')
    neutral_mask = np.isin(char_array, list(neutral_chars))
    mismatches = char_array[:, None, :] != char_array[None, :, :]
    neutrals = neutral_mask[:, None, :] | neutral_mask[None, :, :]
    effective = mismatches & ~neutrals
    return effective.sum(axis=2).astype(float)

def compute_linkage_and_dist(profiles, neutral_chars={'.','N'}, max_tensor_cells=2_50_000_000):
    """Compute distance matrix and hierarchical linkage.

    Returns (headers, dist_matrix, linkage_matrix)."""
    headers = list(profiles.keys())
    n = len(headers)
    dist_matrix = compute_distance_matrix(profiles, neutral_chars=neutral_chars, max_tensor_cells=max_tensor_cells)
    if n < 2:
        return headers, dist_matrix, None
    triu_idx = np.triu_indices(n, k=1)
    distances = dist_matrix[triu_idx]
    Z = linkage(distances, method='average')
    return headers, dist_matrix, Z

def compute_cut_height(Z, max_clusters: int = 10):
    """Heuristic cut height: largest gap in merge distances within cluster count bounds."""
    if Z is None or len(Z) == 0:
        return 0.0, 1
    merge_distances = Z[:, 2]
    gaps = np.diff(merge_distances)
    max_gap = -1.0
    max_idx = -1
    n = Z.shape[0] + 1
    for i, g in enumerate(gaps):
        num_clusters = n - i - 1
        if 2 <= num_clusters <= max_clusters and g > max_gap:
            max_gap = g
            max_idx = i
    if max_idx == -1:
        return float(merge_distances[-1]) if len(merge_distances) else 0.0, 1
    cut_height = merge_distances[max_idx] + (gaps[max_idx] / 2.0)
    num_clusters = n - max_idx - 1
    return float(cut_height), int(num_clusters)

def hierarchical_cluster(profiles, max_clusters=10, neutral_chars={'.','N'}, max_tensor_cells=2_50_000_000):
    """Cluster reads using hierarchical clustering + gap heuristic."""
    if len(profiles) < 2:
        return [set(profiles.keys())], None, None
    headers, dist_matrix, Z = compute_linkage_and_dist(profiles, neutral_chars=neutral_chars, max_tensor_cells=max_tensor_cells)
    if Z is None:
        return [set(headers)], None, dist_matrix
    cut_height, n_clusters = compute_cut_height(Z, max_clusters=max_clusters)
    if n_clusters <= 1:
        log("No significant clustering gap found, using single cluster")
        return [set(headers)], Z, dist_matrix
    log(f"Hierarchical clustering: cutting at height {cut_height:.1f} -> {n_clusters} clusters")
    cluster_labels = fcluster(Z, cut_height, criterion='distance')
    clusters_dict = defaultdict(set)
    for i, label in enumerate(cluster_labels):
        clusters_dict[label].add(headers[i])
    return list(clusters_dict.values()), Z, dist_matrix


def get_parser() -> argparse.ArgumentParser:
    """Return argument parser (facilitates testing/import)."""
    p = argparse.ArgumentParser(description="Detect subclusters from read profiles")
    p.add_argument("profile_file", help="Path to read_profiles.tsv file")
    p.add_argument("outdir", help="Output directory")
    p.add_argument("--viz_out", default=None,
                   help="If set, write visualization to this filename in output directory")
    p.add_argument("--min_cluster_size", type=int, default=5,
                   help="Minimum absolute cluster size")
    p.add_argument("--min_cluster_size_percent", type=float, default=0.0,
                   help="Minimum cluster size as percentage of total reads")
    p.add_argument("--max_clusters", type=int, default=10,
                   help="Maximum number of clusters to detect")
    p.add_argument("--min_variable_positions", type=int, default=3,
                   help="Minimum variable positions required to attempt clustering")
    p.add_argument("--min_reads_to_cluster", type=int, default=None,
                   help="Minimum reads required (default 2 * min_cluster_size)")
    return p

def check_clustering_preconditions(profiles, variable_positions, args, outdir, viz_out_path):
    """Validate variation and read count; handle trivial outputs. Returns True if clustering should proceed."""
    if len(variable_positions) < args.min_variable_positions:
        log(f"Too few variable positions for meaningful clustering (found={len(variable_positions)}, required={args.min_variable_positions})")
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        if viz_out_path:
            write_trivial_profiles_viz(viz_out_path, "No variable positions to visualize")
        return False
    min_reads_required = args.min_reads_to_cluster if args.min_reads_to_cluster is not None else (args.min_cluster_size * 2)
    if len(profiles) < min_reads_required:
        log(f"Too few sequences for clustering ({len(profiles)} < {min_reads_required})")
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        if viz_out_path:
            write_trivial_profiles_viz(viz_out_path, "Not enough reads to cluster")
        return False
    return True

def assign_cluster_labels(valid_clusters):
    """Return dict read_id -> cluster label given list of valid clusters (size-filtered)."""
    cluster_assignments = {}
    if len(valid_clusters) <= 1:
        return cluster_assignments
    cluster_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i, cluster in enumerate(sorted(valid_clusters, key=len, reverse=True)):
        label = cluster_labels[i] if i < len(cluster_labels) else str(i)
        for read_id in sorted(cluster):
            cluster_assignments[read_id] = label
    return cluster_assignments


def main():
    parser = get_parser()
    args = parser.parse_args()
    
    profile_file = Path(args.profile_file)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    viz_out_path = (outdir / args.viz_out) if args.viz_out else None
    
    # Check if profile file exists
    if not profile_file.exists():
        log(f"Error: Profile file not found: {profile_file}")
        sys.exit(1)
    
    # Read profiles
    profiles, variable_positions = read_profiles(profile_file)
    log(f"Loaded {len(profiles)} read profiles with {len(variable_positions)} variable positions")
    
    # Precondition checks (variation & read count)
    if not check_clustering_preconditions(profiles, variable_positions, args, outdir, viz_out_path):
        return
    
    # Cluster reads
    clusters, linkage_matrix, dist_matrix = hierarchical_cluster(profiles, max_clusters=args.max_clusters)
    log(f"Found {len(clusters)} clusters")
    
    # Calculate minimum cluster size
    total_reads = len(profiles)
    min_size_from_percent = int(total_reads * (args.min_cluster_size_percent / 100.0))
    effective_min_size = max(args.min_cluster_size, min_size_from_percent)
    
    log(f"Minimum cluster size: {effective_min_size} reads (absolute={args.min_cluster_size}, percent={args.min_cluster_size_percent}% = {min_size_from_percent} reads)")
    
    # Filter clusters by size
    valid_clusters = [c for c in clusters if len(c) >= effective_min_size]
    log(f"{len(valid_clusters)} clusters meet minimum size threshold")
    
    # Build cluster assignments for visualization (valid clusters only)
    headers = list(profiles.keys())
    cluster_assignments = assign_cluster_labels(valid_clusters)

    # Generate distance heatmap
    distance_heatmap_path = outdir / "distance_heatmap.png"
    try:
        # Fallback to zero matrix if no distances available (e.g., single read)
        if dist_matrix is None:
            dist_matrix = np.zeros((len(headers), len(headers)))
        plot_distance_heatmap(dist_matrix, headers, linkage_matrix, distance_heatmap_path, cluster_assignments=cluster_assignments if cluster_assignments else None)
    except Exception as e:
        log(f"Warning: failed to generate distance heatmap: {e}")
    
    # Produce integrated visualization (with side color bar) if requested
    if viz_out_path:
        try:
            plot_profiles_clustermap(headers, variable_positions, profiles, linkage_matrix, viz_out_path, cluster_assignments=cluster_assignments)
        except Exception as e:
            log(f"Warning: failed to generate integrated viz: {e}")
    
    if len(valid_clusters) < 2:
        log("Only one valid cluster - no subclustering")
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        log("No clustering performed - single cluster")
        return
    
    # Multiple clusters found - write cluster assignments
    log(f"Multiple clusters detected ({len(valid_clusters)} clusters)")
    
    cluster_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    with open(outdir / "cluster_assignments.tsv", "w") as f:
        f.write("read_id\tcluster\n")
        for i, cluster in enumerate(sorted(valid_clusters, key=len, reverse=True)):
            label = cluster_labels[i] if i < len(cluster_labels) else str(i)
            for read_id in sorted(cluster):
                f.write(f"{read_id}\t{label}\n")
    
    log("Wrote cluster assignments to cluster_assignments.tsv")
    
    # Write cluster summary
    with open(outdir / "cluster_summary.tsv", "w") as f:
        f.write("cluster\tsize\n")
        for i, cluster in enumerate(sorted(valid_clusters, key=len, reverse=True)):
            label = cluster_labels[i] if i < len(cluster_labels) else str(i)
            f.write(f"{label}\t{len(cluster)}\n")
    
    # Note: per-request, the primary visualization is the profiles heatmap, already generated above if viz_out set.
    
    log("Cluster detection complete")


if __name__ == "__main__":
    main()
