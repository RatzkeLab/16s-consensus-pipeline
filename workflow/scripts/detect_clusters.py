#!/usr/bin/env python3
"""
Detect if a sample has obvious subclusters.

This script:
1. Identifies variable positions in alignment
2. Creates per-read profiles
3. Performs hierarchical clustering
4. Determines if multiple valid clusters exist
5. Outputs cluster assignments IF clusters are found
"""

import sys
import argparse
from pathlib import Path
from collections import Counter, defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, fcluster


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


def calculate_auto_trim(seqs):
    """Calculate auto-trim values based on longest leading/trailing gaps.
    
    Returns the position of the first non-gap character in the sequence with
    the longest leading gap, and the position of the last non-gap character
    in the sequence with the longest trailing gap.
    
    Args:
        seqs: Dictionary of sequences (aligned)
    
    Returns:
        tuple: (trim_start, trim_end) - number of positions to trim from start and end
    """
    if not seqs:
        return 0, 0
    
    max_leading_gaps = 0
    max_trailing_gaps = 0
    
    for seq in seqs.values():
        # Count leading gaps
        leading = 0
        for char in seq:
            if char == '-':
                leading += 1
            else:
                break
        
        # Count trailing gaps
        trailing = 0
        for char in reversed(seq):
            if char == '-':
                trailing += 1
            else:
                break
        
        max_leading_gaps = max(max_leading_gaps, leading)
        max_trailing_gaps = max(max_trailing_gaps, trailing)
    
    return max_leading_gaps, max_trailing_gaps


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


def hamming_distance(profile1, profile2):
    """Calculate Hamming distance between two profiles."""
    return sum(c1 != c2 for c1, c2 in zip(profile1, profile2))


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
    parser = argparse.ArgumentParser(description="Detect subclusters in sample")
    parser.add_argument("alignment", help="Input alignment FASTA file")
    parser.add_argument("outdir", help="Output directory")
    parser.add_argument("sample", help="Sample name")
    parser.add_argument("--min_agreement", type=float, default=0.8)
    parser.add_argument("--min_cluster_size", type=int, default=5)
    parser.add_argument("--min_cluster_size_percent", type=float, default=0.0)
    parser.add_argument("--max_clusters", type=int, default=10)
    parser.add_argument("--min_variable_positions", type=int, default=3,
                        help="Minimum number of variable positions required to attempt clustering")
    parser.add_argument("--trim_bp", type=int, default=70,
                        help="Number of bp to ignore at start and end of alignment")
    parser.add_argument("--auto_trim", action="store_true",
                        help="Automatically calculate trim values based on longest leading/trailing gaps")
    
    args = parser.parse_args()
    
    # Read alignment
    seqs = read_fasta(args.alignment)
    sys.stderr.write(f"Loaded {len(seqs)} sequences from alignment\n")
    
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    # Check if we have enough sequences
    if len(seqs) < args.min_cluster_size * 2:
        sys.stderr.write(f"Too few sequences for clustering ({len(seqs)} < {args.min_cluster_size * 2})\n")
        # Write empty marker file to indicate no clustering
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        sys.stderr.write("No clustering performed - single cluster\n")
        return
    
    # Identify variable positions
    variable_positions = identify_variable_positions(
        seqs, args.min_agreement, args.trim_bp, args.auto_trim
    )
    sys.stderr.write(f"Found {len(variable_positions)} variable positions\n")
    
    if len(variable_positions) < args.min_variable_positions:
        sys.stderr.write(
            f"Too few variable positions for meaningful clustering "
            f"(found={len(variable_positions)}, required={args.min_variable_positions})\n"
        )
        with open(outdir / "no_clusters.txt", "w") as f:
            f.write("single_cluster\n")
        sys.stderr.write("No clustering performed - single cluster\n")
        return
    
    # Create profiles
    profiles = create_read_profiles(seqs, variable_positions)
    
    # Cluster reads
    clusters, linkage_matrix, dist_matrix = hierarchical_cluster(profiles, max_clusters=args.max_clusters)
    sys.stderr.write(f"Found {len(clusters)} clusters\n")
    
    # Calculate minimum cluster size
    total_reads = len(seqs)
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
