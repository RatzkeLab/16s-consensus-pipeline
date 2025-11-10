#!/usr/bin/env python3
"""
Build a symmetric edit distance matrix from pairwise TSV, cluster it, and plot a heatmap.

Inputs (from snakemake):
- snakemake.input.distances: TSV with columns [Sequence1, Sequence2, EditDistance]

Outputs (from snakemake):
- snakemake.output.matrix: TSV square matrix with clustered row/column order
- snakemake.output.heatmap: PNG/PDF heatmap image of the clustered matrix

Clustering: SciPy hierarchical clustering (average linkage) on the condensed distance matrix.
"""

import sys
from pathlib import Path
from typing import Union, List

import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform


def load_pairwise_tsv(path: Union[str, Path]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    # Basic validation
    expected_cols = {"Sequence1", "Sequence2", "EditDistance"}
    if not expected_cols.issubset(df.columns):
        raise ValueError(f"pairwise TSV missing columns. Found: {df.columns.tolist()}, expected at least: {sorted(expected_cols)}")
    return df


def build_symmetric_matrix(df: pd.DataFrame) -> pd.DataFrame:
    # Get unique labels
    labels = sorted(set(df["Sequence1"]).union(df["Sequence2"]))
    n = len(labels)
    mat = pd.DataFrame(np.zeros((n, n), dtype=float), index=labels, columns=labels)

    # Fill from upper-tri list
    for _, row in df.iterrows():
        a = row["Sequence1"]
        b = row["Sequence2"]
        d = float(row["EditDistance"]) if pd.notna(row["EditDistance"]) else np.nan
        mat.at[a, b] = d
        mat.at[b, a] = d

    # Ensure diagonal is zero
    for i in labels:
        mat.at[i, i] = 0.0

    # If any NaNs, fill with max observed distance (conservative) to keep clustering robust
    if mat.isna().any().any():
        max_d = np.nanmax(mat.values)
        if np.isnan(max_d):
            max_d = 0.0
        mat = mat.fillna(max_d)

    return mat


def cluster_order(distance_matrix: pd.DataFrame) -> tuple:
    """Return clustered order and linkage matrix."""
    n = distance_matrix.shape[0]
    if n <= 2:
        # No meaningful clustering for <=2; keep original order
        return list(distance_matrix.index), None

    # Convert to condensed form required by linkage
    condensed = squareform(distance_matrix.values, checks=False)
    Z = linkage(condensed, method="average")
    return None, Z  # seaborn will handle ordering


def plot_heatmap(dm: pd.DataFrame, out_path: Union[str, Path], linkage_matrix=None) -> None:
    """Plot heatmap with dendrogram using seaborn clustermap and overlaid distance values."""
    n = dm.shape[0]
    
    if linkage_matrix is not None and n > 2:
        # Use seaborn clustermap with precomputed linkage
        # Calculate figure size based on number of sequences
        figsize = (max(10, 0.5 * n + 3), max(8, 0.5 * n + 1))
        
        g = sns.clustermap(dm, 
                          row_linkage=linkage_matrix,
                          col_linkage=linkage_matrix,
                          cmap="viridis",
                          figsize=figsize,
                          cbar_kws={'label': 'Edit distance'},
                          dendrogram_ratio=0.15,
                          xticklabels=True,
                          yticklabels=True)
        
        # Adjust tick label sizes
        g.ax_heatmap.tick_params(axis='both', labelsize=max(6, min(10, 80 / n if n > 0 else 10)))
        plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
        plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
        
        # Overlay text values on heatmap - bold and white
        vmin, vmax = dm.values.min(), dm.values.max()
        
        # Get the reordered indices from the clustermap
        row_order = g.dendrogram_row.reordered_ind
        col_order = g.dendrogram_col.reordered_ind
        
        for i in range(n):
            for j in range(n):
                # Get value from reordered matrix
                val = dm.values[row_order[i], col_order[j]]
                # Format: show integer if whole number, else one decimal
                text = f"{int(val)}" if val == int(val) else f"{val:.1f}"
                g.ax_heatmap.text(j + 0.5, i + 0.5, text, 
                                 ha="center", va="center",
                                 color="white", 
                                 fontsize=max(5, min(9, 70 / n if n > 0 else 9)),
                                 weight='bold')
        
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close()
    else:
        # Simple heatmap without clustering for small matrices
        fig, ax = plt.subplots(figsize=(max(8, 0.5 * n + 2), max(6, 0.5 * n + 1)))
        im = ax.imshow(dm.values, cmap="viridis", aspect="auto")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Edit distance", rotation=270, labelpad=20)
        
        ax.set_xticks(np.arange(n))
        ax.set_yticks(np.arange(n))
        ax.set_xticklabels(dm.columns, rotation=90, fontsize=max(6, min(10, 80 / n if n > 0 else 10)))
        ax.set_yticklabels(dm.index, fontsize=max(6, min(10, 80 / n if n > 0 else 10)))
        
        # Overlay text - bold and white
        for i in range(n):
            for j in range(n):
                val = dm.values[i, j]
                text = f"{int(val)}" if val == int(val) else f"{val:.1f}"
                ax.text(j, i, text, ha="center", va="center",
                       color="white", fontsize=max(5, min(9, 70 / n if n > 0 else 9)),
                       weight='bold')
        
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close()


def main():
    distances_tsv = snakemake.input.distances  # type: ignore[name-defined]
    matrix_tsv = snakemake.output.matrix       # type: ignore[name-defined]
    heatmap_img = snakemake.output.heatmap     # type: ignore[name-defined]

    df = load_pairwise_tsv(distances_tsv)
    mat = build_symmetric_matrix(df)

    # Cluster and get linkage matrix
    _, linkage_matrix = cluster_order(mat)

    # Write matrix (original order)
    out_p = Path(matrix_tsv)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(out_p, sep="\t")

    # Plot heatmap (seaborn will reorder)
    heatmap_p = Path(heatmap_img)
    heatmap_p.parent.mkdir(parents=True, exist_ok=True)
    plot_heatmap(mat, heatmap_p, linkage_matrix)

    print(f"Wrote clustered distance matrix to {out_p}", file=sys.stderr)
    print(f"Wrote heatmap to {heatmap_p}", file=sys.stderr)


if __name__ == "__main__":
    main()
