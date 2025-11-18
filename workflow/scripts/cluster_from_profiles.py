#!/usr/bin/env python3
"""
Detect if a sample has obvious subclusters based on pre-generated read profiles.

This script:
1. Reads per-read profiles from upstream step
2. Performs hierarchical clustering
3. Determines if multiple valid clusters exist
4. Outputs cluster assignments IF clusters are found
"""

import sys
import argparse
from pathlib import Path
from collections import defaultdict
import numpy as np
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
    """Read read profiles from TSV file using vectorized numpy operations.
    
    Returns:
        Tuple of (headers, profiles_array, variable_positions)
        - headers: list of read IDs in file order
        - profiles_array: (n_reads, n_positions) numpy array of characters (dtype='S1')
        - variable_positions: list of position numbers (parsed from header)
    """
    with open(profile_file) as f:
        # Parse header to get variable positions
        header = f.readline().strip().split('\t')
        variable_positions = [int(pos.replace('pos_', '')) for pos in header[1:]]
        
        # Read all lines at once and filter empty
        lines = [line.strip() for line in f if line.strip()]
    
    if not lines:
        return [], np.zeros((0, len(variable_positions)), dtype='S1'), variable_positions
    
    # Vectorized split: use numpy's vectorized string operations
    # Split all lines in batch
    split_lines = [line.split('\t') for line in lines]
    
    # Convert to numpy array for vectorized indexing
    data = np.array(split_lines, dtype=str)
    headers = data[:, 0].tolist()
    profiles_array = data[:, 1:].astype('S1')  # ASCII bytes for memory efficiency
    
    return headers, profiles_array, variable_positions


def validate_profiles(headers, profiles_array, variable_positions):
    """Validate profile consistency and completeness.
    
    Args:
        headers: list of read IDs
        profiles_array: (n_reads, n_positions) numpy array
        variable_positions: list of expected position numbers
    
    Raises:
        ValueError: if profiles have inconsistent lengths or don't match variable_positions
    
    Returns:
        None (raises on error)
    """
    expected_length = len(variable_positions)
    n_reads = len(headers)
    
    if expected_length == 0:
        log("Warning: No variable positions in profile file")
        return
    
    if n_reads == 0:
        log("Warning: No profiles found in file")
        return
    
    # Vectorized validation: check array shape
    if profiles_array.shape != (n_reads, expected_length):
        actual_cols = profiles_array.shape[1] if len(profiles_array.shape) > 1 else 0
        error_msg = f"Profile shape mismatch: expected ({n_reads}, {expected_length}), got {profiles_array.shape}\n"
        error_msg += f"All profiles must have exactly {expected_length} positions.\n"
        raise ValueError(error_msg)
    
    log(f"Profile validation passed: {n_reads} reads × {expected_length} positions")


def condensed_distances(char_array, neutral_mask):
    """Compute condensed pairwise distance vector (upper triangle, 1D).
    
    Returns 1D array of length n*(n-1)/2 suitable for scipy.linkage.
    Distance ignores neutral positions and counts mismatches.
    """
    n, m = char_array.shape
    if n < 2:
        return np.array([])
    # Vectorized pairwise comparison for upper triangle
    # Build index pairs for condensed form
    i_idx, j_idx = np.triu_indices(n, k=1)
    # Broadcasting: compare all pairs
    mismatches = char_array[i_idx, :] != char_array[j_idx, :]  # (num_pairs, m)
    neutrals = neutral_mask[i_idx, :] | neutral_mask[j_idx, :]  # (num_pairs, m)
    effective = mismatches & ~neutrals
    distances = effective.sum(axis=1).astype(float)
    return distances


def square_from_condensed(condensed, n):
    """Build square distance matrix from condensed 1D vector."""
    if n < 2:
        return np.zeros((n, n))
    dist_matrix = np.zeros((n, n), dtype=float)
    i_idx, j_idx = np.triu_indices(n, k=1)
    dist_matrix[i_idx, j_idx] = condensed
    dist_matrix[j_idx, i_idx] = condensed  # Symmetric
    return dist_matrix


def compute_linkage_and_dist(headers, profiles_array, neutral_chars={'.','N'}, linkage_method='average'):
    """Compute condensed distances and hierarchical linkage.
    
    Args:
        headers: list of read IDs
        profiles_array: (n_reads, n_positions) numpy array (dtype='S1')
        neutral_chars: set/iterable of characters to ignore in distance
        linkage_method: scipy linkage method ('average', 'single', 'complete', 'ward', etc.)
    
    Returns (condensed_distances, linkage_matrix).
    Condensed distances are 1D (n*(n-1)/2); use square_from_condensed for plotting.
    """
    n = len(headers)
    if n < 2:
        return np.array([]), None
    
    # Build neutral mask directly from profiles_array
    neutral_list = [s.encode('ascii') if isinstance(s, str) else s for s in neutral_chars]
    neutral_mask = np.isin(profiles_array, neutral_list)
    
    distances = condensed_distances(profiles_array, neutral_mask)
    Z = linkage(distances, method=linkage_method)
    return distances, Z

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

def hierarchical_cluster(headers, profiles_array, max_clusters=10, neutral_chars={'.','N'}, linkage_method='average'):
    """Cluster reads using hierarchical clustering + gap heuristic.
    
    Args:
        headers: list of read IDs
        profiles_array: (n_reads, n_positions) numpy array (dtype='S1')
        max_clusters: maximum number of clusters to detect
        neutral_chars: set/iterable of characters to ignore in distance (configurable via CLI)
        linkage_method: scipy linkage method (configurable via CLI)
    
    Returns (clusters, linkage_matrix, condensed_distances):
        - clusters: list of sets of read IDs
        - linkage_matrix: scipy linkage matrix (or None)
        - condensed_distances: 1D distance vector (or empty array)
    """
    if len(headers) < 2:
        return [set(headers)], None, np.array([])
    condensed_dist, Z = compute_linkage_and_dist(headers, profiles_array, neutral_chars=neutral_chars, linkage_method=linkage_method)
    if Z is None:
        return [set(headers)], None, condensed_dist
    cut_height, n_clusters = compute_cut_height(Z, max_clusters=max_clusters)
    if n_clusters <= 1:
        log("No significant clustering gap found, using single cluster")
        return [set(headers)], Z, condensed_dist
    log(f"Hierarchical clustering: cutting at height {cut_height:.1f} -> {n_clusters} clusters")
    cluster_labels = fcluster(Z, cut_height, criterion='distance')
    clusters_dict = defaultdict(set)
    for i, label in enumerate(cluster_labels):
        clusters_dict[label].add(headers[i])
    return list(clusters_dict.values()), Z, condensed_dist


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
    p.add_argument("--neutral_chars", type=str, default=".N",
                   help="Characters to treat as neutral (ignored in distance). Default: '.N'")
    p.add_argument("--linkage_method", type=str, default="average",
                   choices=["average", "single", "complete", "weighted", "centroid", "median", "ward"],
                   help="Linkage method for hierarchical clustering. Default: 'average'")
    return p

def check_clustering_preconditions(n_reads, variable_positions, args, outdir, viz_out_path):
    """Validate variation and read count; handle trivial outputs. Returns True if clustering should proceed."""
    if len(variable_positions) < args.min_variable_positions:
        log(f"Too few variable positions for meaningful clustering (found={len(variable_positions)}, required={args.min_variable_positions})")
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        if viz_out_path:
            write_trivial_profiles_viz(viz_out_path, "No variable positions to visualize")
        return False
    min_reads_required = args.min_reads_to_cluster if args.min_reads_to_cluster is not None else (args.min_cluster_size * 2)
    if n_reads < min_reads_required:
        log(f"Too few sequences for clustering ({n_reads} < {min_reads_required})")
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
    
    # Parse neutral characters (convert string to set)
    neutral_chars = set(args.neutral_chars)
    log(f"Neutral characters (ignored in distance): {sorted(neutral_chars)}")
    log(f"Linkage method: {args.linkage_method}")
    
    # Check if profile file exists
    if not profile_file.exists():
        log(f"Error: Profile file not found: {profile_file}")
        sys.exit(1)
    
    # Read profiles
    headers, profiles_array, variable_positions = read_profiles(profile_file)
    log(f"Loaded {len(headers)} read profiles with {len(variable_positions)} variable positions")
    
    # Validate profile consistency
    try:
        validate_profiles(headers, profiles_array, variable_positions)
    except ValueError as e:
        log(f"Error: Profile validation failed")
        log(str(e))
        sys.exit(1)
    
    # Precondition checks (variation & read count)
    if not check_clustering_preconditions(len(headers), variable_positions, args, outdir, viz_out_path):
        return
    
    # Cluster reads
    clusters, linkage_matrix, condensed_dist = hierarchical_cluster(
        headers,
        profiles_array,
        max_clusters=args.max_clusters, 
        neutral_chars=neutral_chars,
        linkage_method=args.linkage_method
    )
    log(f"Found {len(clusters)} clusters")
    
    # Calculate minimum cluster size
    total_reads = len(headers)
    min_size_from_percent = int(total_reads * (args.min_cluster_size_percent / 100.0))
    effective_min_size = max(args.min_cluster_size, min_size_from_percent)
    
    log(f"Minimum cluster size: {effective_min_size} reads (absolute={args.min_cluster_size}, percent={args.min_cluster_size_percent}% = {min_size_from_percent} reads)")
    
    # Filter clusters by size
    valid_clusters = [c for c in clusters if len(c) >= effective_min_size]
    log(f"{len(valid_clusters)} clusters meet minimum size threshold")
    
    # Build cluster assignments for visualization (valid clusters only)
    cluster_assignments = assign_cluster_labels(valid_clusters)

    # Generate distance heatmap (build square matrix only when needed for plotting)
    distance_heatmap_path = outdir / "distance_heatmap.png"
    try:
        n = len(headers)
        if len(condensed_dist) == 0:
            dist_matrix = np.zeros((n, n))
        else:
            dist_matrix = square_from_condensed(condensed_dist, n)
        plot_distance_heatmap(dist_matrix, headers, linkage_matrix, distance_heatmap_path, cluster_assignments=cluster_assignments if cluster_assignments else None)
    except Exception as e:
        log(f"Warning: failed to generate distance heatmap: {e}")
    
    # Produce integrated visualization (with side color bar) if requested
    if viz_out_path:
        try:
            plot_profiles_clustermap(headers, profiles_array, variable_positions, linkage_matrix, viz_out_path, cluster_assignments=cluster_assignments)
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
