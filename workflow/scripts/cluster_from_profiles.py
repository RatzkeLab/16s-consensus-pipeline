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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster


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

    Rules:
    - 'N' is treated as neutral (wildcard) and does not contribute to distance
        when present on either side.
    - All other characters are compared directly; mismatch counts as 1.
    - Gaps ('-') are penalized only for the first position of a gap run when
        used together with compress_gap_runs() which marks continuation gaps as 'N'.
    """
    dist = 0
    for c1, c2 in zip(profile1, profile2):
        # Ignore positions where either side is neutral
        if c1 == 'N' or c2 == 'N':
            continue
        if c1 != c2:
            dist += 1
    return dist


def hierarchical_cluster(profiles, max_clusters=10):
    """
    Cluster reads using hierarchical clustering with automatic cut height detection.
    
    Returns:
        Tuple of (clusters, linkage_matrix, distance_matrix)
        - clusters: List of clusters, where each cluster is a set of read headers
        - linkage_matrix: scipy linkage matrix for dendrogram plotting
        - distance_matrix: full distance matrix (square) as numpy array
    """
    if len(profiles) < 2:
        return [set(profiles.keys())], None, None
    
    headers = list(profiles.keys())
    n = len(headers)
    
    # Build distance matrix (condensed form for linkage)
    distances = []
    for i in range(n):
        for j in range(i + 1, n):
            dist = hamming_distance(profiles[headers[i]], profiles[headers[j]])
            distances.append(dist)
    
    # Build full square distance matrix for visualization
    dist_matrix = np.zeros((n, n))
    idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            dist_matrix[i, j] = distances[idx]
            dist_matrix[j, i] = distances[idx]
            idx += 1
    
    # Perform hierarchical clustering using average linkage
    Z = linkage(distances, method='average')
    
    # Find optimal cut height by looking for largest gap in merge distances
    merge_distances = Z[:, 2]
    gaps = np.diff(merge_distances)
    
    # Find the largest gap
    max_gap_idx = -1
    max_gap = -1
    
    for idx in range(len(gaps)):
        num_clusters = n - idx - 1
        if 2 <= num_clusters <= max_clusters and gaps[idx] > max_gap:
            max_gap = gaps[idx]
            max_gap_idx = idx
    
    if max_gap_idx == -1:
        sys.stderr.write("No significant clustering gap found, using single cluster\n")
        return [set(headers)], Z, dist_matrix
    
    # Cut dendrogram at the merge just before the largest gap
    cut_height = merge_distances[max_gap_idx] + (gaps[max_gap_idx] / 2)
    num_clusters = n - max_gap_idx - 1
    
    sys.stderr.write(f"Hierarchical clustering: cutting at height {cut_height:.1f} "
                    f"(gap={max_gap:.1f}) -> {num_clusters} clusters\n")
    
    # Get cluster labels
    cluster_labels = fcluster(Z, cut_height, criterion='distance')
    
    # Group headers by cluster
    clusters_dict = defaultdict(set)
    for i, label in enumerate(cluster_labels):
        clusters_dict[label].add(headers[i])
    
    return list(clusters_dict.values()), Z, dist_matrix


def plot_cluster_heatmap(dist_matrix, headers, linkage_matrix, outpath, cluster_assignments=None):
    """
    Plot a clustermap of the distance matrix.
    
    Args:
        dist_matrix: Square distance matrix (numpy array)
        headers: List of read IDs corresponding to matrix rows/columns
        linkage_matrix: scipy linkage matrix for dendrograms
        outpath: Output file path for the heatmap
        cluster_assignments: Optional dict mapping read_id -> cluster_label for color annotation
    """
    import pandas as pd
    
    # Create dataframe from distance matrix
    df = pd.DataFrame(dist_matrix, index=headers, columns=headers)
    
    # Determine figure size based on number of reads
    n = len(headers)
    figsize = (max(10, 0.15 * n), max(8, 0.15 * n))
    
    # Create row colors if cluster assignments provided
    row_colors = None
    if cluster_assignments:
        # Map cluster labels to colors
        # Only include reads that have cluster assignments
        unique_clusters = sorted(set(cluster_assignments.values()))
        color_palette = sns.color_palette("tab10", n_colors=len(unique_clusters))
        cluster_color_map = {cluster: color_palette[i] for i, cluster in enumerate(unique_clusters)}
        
        # Default color for unassigned reads (gray)
        default_color = (0.7, 0.7, 0.7)
        row_colors = [cluster_color_map.get(cluster_assignments.get(h), default_color) for h in headers]
    
    # Create clustermap
    try:
        g = sns.clustermap(
            df,
            row_linkage=linkage_matrix,
            col_linkage=linkage_matrix,
            cmap="viridis",
            figsize=figsize,
            cbar_kws={'label': 'Hamming distance'},
            dendrogram_ratio=0.15,
            xticklabels=False,  # Too many labels for reads
            yticklabels=False,
            row_colors=row_colors,
            colors_ratio=0.03 if row_colors else None
        )
        
        g.ax_heatmap.set_xlabel("Reads")
        g.ax_heatmap.set_ylabel("Reads")
        
        plt.savefig(outpath, dpi=150, bbox_inches='tight')
        plt.close()
        sys.stderr.write(f"Wrote cluster heatmap to {outpath}\n")
    except Exception as e:
        sys.stderr.write(f"Warning: Could not create clustermap: {e}\n")


def main():
    parser = argparse.ArgumentParser(description="Detect subclusters from read profiles")
    parser.add_argument("profile_file", help="Path to read_profiles.tsv file")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("--min_cluster_size", type=int, default=5,
                        help="Minimum absolute cluster size")
    parser.add_argument("--min_cluster_size_percent", type=float, default=0.0,
                        help="Minimum cluster size as percentage of total reads")
    parser.add_argument("--max_clusters", type=int, default=10,
                        help="Maximum number of clusters to detect")
    parser.add_argument("--min_variable_positions", type=int, default=3,
                        help="Minimum number of variable positions required to attempt clustering")
    
    args = parser.parse_args()
    
    profile_file = Path(args.profile_file)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Check if profile file exists
    if not profile_file.exists():
        sys.stderr.write(f"Error: Profile file not found: {profile_file}\n")
        sys.exit(1)
    
    # Read profiles
    profiles, variable_positions = read_profiles(profile_file)
    sys.stderr.write(f"Loaded {len(profiles)} read profiles with {len(variable_positions)} variable positions\n")
    
    # Check if we have enough variable positions for clustering
    if len(variable_positions) < args.min_variable_positions:
        sys.stderr.write(
            f"Too few variable positions for meaningful clustering "
            f"(found={len(variable_positions)}, required={args.min_variable_positions})\n"
        )
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        sys.stderr.write("No clustering performed - insufficient variation\n")
        return
    
    # Check if we have enough sequences for clustering
    if len(profiles) < args.min_cluster_size * 2:
        sys.stderr.write(f"Too few sequences for clustering ({len(profiles)} < {args.min_cluster_size * 2})\n")
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        sys.stderr.write("No clustering performed - single cluster\n")
        return
    
    # Cluster reads
    clusters, linkage_matrix, dist_matrix = hierarchical_cluster(profiles, max_clusters=args.max_clusters)
    sys.stderr.write(f"Found {len(clusters)} clusters\n")
    
    # Calculate minimum cluster size
    total_reads = len(profiles)
    min_size_from_percent = int(total_reads * (args.min_cluster_size_percent / 100.0))
    effective_min_size = max(args.min_cluster_size, min_size_from_percent)
    
    sys.stderr.write(f"Minimum cluster size: {effective_min_size} reads "
                    f"(absolute={args.min_cluster_size}, "
                    f"percent={args.min_cluster_size_percent}% = {min_size_from_percent} reads)\n")
    
    # Filter clusters by size
    valid_clusters = [c for c in clusters if len(c) >= effective_min_size]
    sys.stderr.write(f"{len(valid_clusters)} clusters meet minimum size threshold\n")
    
    if len(valid_clusters) < 2:
        sys.stderr.write("Only one valid cluster - no subclustering\n")
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        sys.stderr.write("No clustering performed - single cluster\n")
        
        # Still plot heatmap for QC (even if no clusters)
        if linkage_matrix is not None and dist_matrix is not None:
            headers = list(profiles.keys())
            plot_cluster_heatmap(
                dist_matrix, 
                headers, 
                linkage_matrix, 
                outdir / "cluster_heatmap.png"
            )
        return
    
    # Multiple clusters found - write cluster assignments
    sys.stderr.write(f"Multiple clusters detected ({len(valid_clusters)} clusters)\n")
    
    cluster_labels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    cluster_assignments = {}
    
    with open(outdir / "cluster_assignments.tsv", "w") as f:
        f.write("read_id\tcluster\n")
        for i, cluster in enumerate(sorted(valid_clusters, key=len, reverse=True)):
            label = cluster_labels[i] if i < len(cluster_labels) else str(i)
            for read_id in sorted(cluster):
                f.write(f"{read_id}\t{label}\n")
                cluster_assignments[read_id] = label
    
    sys.stderr.write(f"Wrote cluster assignments to cluster_assignments.tsv\n")
    
    # Write cluster summary
    with open(outdir / "cluster_summary.tsv", "w") as f:
        f.write("cluster\tsize\n")
        for i, cluster in enumerate(sorted(valid_clusters, key=len, reverse=True)):
            label = cluster_labels[i] if i < len(cluster_labels) else str(i)
            f.write(f"{label}\t{len(cluster)}\n")
    
    # Plot cluster heatmap with cluster annotations
    if linkage_matrix is not None and dist_matrix is not None:
        headers = list(profiles.keys())
        plot_cluster_heatmap(
            dist_matrix,
            headers,
            linkage_matrix,
            outdir / "cluster_heatmap.png",
            cluster_assignments
        )
    
    sys.stderr.write("Cluster detection complete\n")


if __name__ == "__main__":
    main()
