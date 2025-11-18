#!/usr/bin/env python3
"""Visualization helpers for clustering read profiles.

Separated from cluster_from_profiles.py to reduce bulk and isolate plotting concerns.
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


def log(msg: str) -> None:
    """Local logger (stderr). Keep lightweight to avoid circular imports."""
    sys.stderr.write(f"{msg}\n")

# -------------------- Distance heatmap --------------------

def plot_distance_heatmap(dist_matrix, headers, linkage_matrix, outpath, cluster_assignments=None):
    """Generate a heatmap of pairwise distances between reads.

    Args:
        dist_matrix: Square distance matrix (numpy array)
        headers: List of read IDs corresponding to matrix rows/columns
        linkage_matrix: scipy linkage matrix for dendrograms
        outpath: Output file path for the heatmap
        cluster_assignments: Optional dict mapping read_id -> cluster_label for color annotation
    """
    import pandas as pd

    if len(headers) == 0:
        log("No reads to visualize in distance heatmap")
        return

    df = pd.DataFrame(dist_matrix, index=headers, columns=headers)
    n = len(headers)
    figsize = (max(10, 0.15 * n), max(8, 0.15 * n))

    row_colors = None
    if cluster_assignments:
        unique_clusters = sorted(set(cluster_assignments.values()))
        palette_name = "tab10" if len(unique_clusters) <= 10 else "tab20"
        color_palette = sns.color_palette(palette_name, n_colors=len(unique_clusters))
        cluster_color_map = {cluster: color_palette[i] for i, cluster in enumerate(unique_clusters)}
        default_color = (0.7, 0.7, 0.7)
        row_colors = [cluster_color_map.get(cluster_assignments.get(h), default_color) for h in headers]

    labels_trunc = [h[:15] for h in headers]

    try:
        g = sns.clustermap(
            df,
            row_linkage=linkage_matrix,
            col_linkage=linkage_matrix,
            cmap="viridis",
            figsize=figsize,
            cbar_kws={'label': 'Hamming distance'},
            dendrogram_ratio=0.15,
            xticklabels=False,
            yticklabels=True,
            row_colors=row_colors,
            colors_ratio=0.03 if row_colors else None
        )
        g.ax_heatmap.set_xlabel("Reads")
        g.ax_heatmap.set_ylabel("Reads")
        g.ax_heatmap.set_title("Pairwise distance heatmap")
        try:
            if hasattr(g, 'dendrogram_row') and hasattr(g.dendrogram_row, 'reordered_ind'):
                order = g.dendrogram_row.reordered_ind
                ordered_labels = [labels_trunc[i] for i in order]
                g.ax_heatmap.set_yticklabels(ordered_labels, rotation=0)
            else:
                g.ax_heatmap.set_yticklabels(labels_trunc, rotation=0)
            g.ax_heatmap.tick_params(axis='y', labelsize=max(4, min(9, int(120 / max(n, 1)))))
        except Exception as e:
            log(f"Warning: failed to set y-axis labels: {e}")
        plt.savefig(outpath, dpi=150, bbox_inches='tight')
        plt.close()
        log(f"Wrote distance heatmap to {outpath}")
    except Exception as e:
        log(f"Warning: Could not create distance heatmap: {e}")

# -------------------- Profile heatmap --------------------

def profiles_to_numeric_matrix(profiles_array):
    """Convert categorical profiles array to a numeric matrix and return a colormap/norm.

    Args:
        profiles_array: (n_reads, n_positions) numpy array of ASCII bytes (dtype='S1')

    Returns (matrix (n_reads x n_positions), cmap, norm, legend_labels)
    """
    from matplotlib.colors import ListedColormap, BoundaryNorm
    symbols = ['A', 'C', 'G', 'T', '-', '.', 'N']
    sym_bytes = [s.encode('ascii') for s in symbols]
    sym_index = {s: i for i, s in enumerate(sym_bytes)}
    unknown_idx = len(symbols)
    
    # Vectorized mapping
    mat = np.full(profiles_array.shape, unknown_idx, dtype=int)
    for i, byte_sym in enumerate(sym_bytes):
        mat[profiles_array == byte_sym] = i
    
    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#bdbdbd", "#aaaaaa", "#000000"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(-0.5, len(colors) + 0.5), len(colors))
    return mat, cmap, norm, symbols + ['?']


def write_trivial_profiles_viz(out_png, message):
    """Write a simple explanatory figure when no data to visualize."""
    try:
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, message, ha='center', va='center')
        ax.axis('off')
        out_path = Path(out_png)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        log(f"Wrote trivial profiles viz to {out_path}")
    except Exception as e:
        log(f"Warning: failed to write trivial profiles viz: {e}")


def overlay_profile_text(g, profiles_array, headers):
    """Overlay per-position characters on existing clustermap heatmap.
    
    Args:
        g: seaborn clustermap object
        profiles_array: (n_reads, n_positions) numpy array (dtype='S1')
        headers: list of read IDs
    """
    try:
        h = len(headers)
        if h == 0:
            return
        row_order = g.dendrogram_row.reordered_ind if hasattr(g, 'dendrogram_row') else list(range(h))
        w = profiles_array.shape[1] if profiles_array.size > 0 else 0
        fontsize = max(3, min(8, 50 / max(h, w, 1)))
        for i in range(h):
            orig_row_idx = row_order[i] if i < len(row_order) else i
            if orig_row_idx >= len(headers):
                continue
            prof_bytes = profiles_array[orig_row_idx, :]
            for j, byte_letter in enumerate(prof_bytes):
                letter = byte_letter.decode('ascii') if isinstance(byte_letter, bytes) else str(byte_letter)
                g.ax_heatmap.text(j + 0.5, i + 0.5, letter,
                                  ha="center", va="center",
                                  color="white" if letter not in ['.', 'N'] else "black",
                                  fontsize=fontsize,
                                  weight='bold')
    except Exception as e:
        log(f"Warning: overlay text failed: {e}")


def add_cluster_legend(g, cluster_assignments, label_to_color):
    """Attach cluster legend to figure if assignments present."""
    if not cluster_assignments or not label_to_color:
        return
    try:
        from matplotlib.patches import Patch
        labels_present = sorted(set(cluster_assignments.values()))
        handles = [Patch(facecolor=label_to_color[lab], edgecolor='none', label=str(lab))
                   for lab in labels_present if lab in label_to_color]
        if handles:
            g.fig.legend(handles, [h.get_label() for h in handles], loc='upper left',
                         title='Cluster', frameon=False)
    except Exception as e:
        log(f"Warning: could not add cluster legend: {e}")


def add_symbol_legend(g, legend_labels):
    """Attach symbol legend (profile character color mapping)."""
    if not legend_labels:
        return
    try:
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=g.cmap(i), edgecolor='none', label=legend_labels[i])
                   for i in range(len(legend_labels))]
        g.fig.legend(handles, legend_labels, loc='upper right', title='Symbol', frameon=False)
    except Exception:
        pass


def plot_profiles_clustermap(headers, profiles_array, variable_positions, linkage_matrix, out_png, cluster_assignments=None):
    """Plot profile categorical clustermap with dendrogram and optional cluster side colors.
    
    Args:
        headers: list of read IDs
        profiles_array: (n_reads, n_positions) numpy array (dtype='S1')
        variable_positions: list of position numbers
        linkage_matrix: scipy linkage matrix
        out_png: output file path
        cluster_assignments: optional dict mapping read_id -> cluster_label
    """
    from matplotlib.patches import Patch  # noqa: F401 (used indirectly)
    if len(headers) == 0 or len(variable_positions) == 0:
        write_trivial_profiles_viz(out_png, "No reads to visualize" if len(headers) == 0 else "No variable positions to visualize")
        return
    try:
        mat, cmap, norm, legend_labels = profiles_to_numeric_matrix(profiles_array)
    except Exception as e:
        log(f"Warning: failed to build numeric matrix: {e}")
        return
    import pandas as pd
    df = pd.DataFrame(mat, index=headers, columns=[f"pos_{p}" for p in variable_positions])
    h, w = df.shape
    if h == 0 or w == 0:
        write_trivial_profiles_viz(out_png, "Empty profiles frame")
        return
    figsize = (max(6, min(18, 0.12 * w + 4)), max(6, min(18, 0.12 * h + 4)))
    row_colors = None
    label_to_color = {}
    if cluster_assignments:
        try:
            labels_present = sorted(set(cluster_assignments.values()))
            palette_name = "tab20" if len(labels_present) <= 20 else "hsv"
            pal = sns.color_palette(palette_name, n_colors=len(labels_present))
            label_to_color = {lab: pal[i] for i, lab in enumerate(labels_present)}
            default_color = (0.7, 0.7, 0.7)
            row_colors = [label_to_color.get(cluster_assignments.get(h), default_color) for h in headers]
        except Exception as e:
            log(f"Warning: failed to build cluster colors: {e}")
            row_colors = None
    try:
        g = sns.clustermap(
            df,
            row_linkage=linkage_matrix if linkage_matrix is not None else None,
            col_cluster=False,
            cmap=cmap,
            norm=norm,
            xticklabels=True,
            yticklabels=True,
            figsize=figsize,
            dendrogram_ratio=0.15,
            cbar_pos=None,
            row_colors=row_colors,
            colors_ratio=0.015 if row_colors is not None else None,
        )
    except Exception as e:
        log(f"Warning: clustermap failed ({e}); falling back to imshow heatmap")
        fig, ax = plt.subplots(figsize=figsize)
        ax.imshow(mat, aspect='auto', cmap=cmap, norm=norm, interpolation='nearest')
        ax.set_xlabel('Variable positions')
        ax.set_ylabel('Reads')
        try:
            labels_trunc = [h[:15] for h in headers]
            ax.set_yticks(range(len(labels_trunc)))
            ax.set_yticklabels(labels_trunc, fontsize=max(4, min(9, int(120 / max(len(labels_trunc), 1)))))
        except Exception:
            pass
        ax.set_title("Profiles heatmap")
        out_path = Path(out_png)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close(fig)
        log(f"Wrote fallback heatmap to {out_path}")
        return
    # Axis adjustments
    try:
        g.ax_heatmap.tick_params(axis='x', labelsize=max(4, min(8, 60 / max(w, 1))), rotation=90)
    except Exception as e:
        log(f"Warning: failed to adjust x-axis labels: {e}")
    try:
        labels_trunc = [h[:15] for h in headers]
        row_order = g.dendrogram_row.reordered_ind if (linkage_matrix is not None and hasattr(g, 'dendrogram_row')) else list(range(h))
        ordered_labels = [labels_trunc[i] for i in row_order]
        g.ax_heatmap.set_yticklabels(ordered_labels, rotation=0)
        g.ax_heatmap.tick_params(axis='y', labelsize=max(4, min(9, int(120 / max(h, 1)))))
    except Exception as e:
        log(f"Warning: failed to set y-axis labels: {e}")
    overlay_profile_text(g, profiles_array, headers)
    try:
        g.ax_heatmap.set_title("Profiles heatmap with dendrogram", pad=12)
    except Exception:
        pass
    add_cluster_legend(g, cluster_assignments, label_to_color)
    add_symbol_legend(g, legend_labels)
    out_path = Path(out_png)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close()
        log(f"Wrote profile clustermap to {out_path}")
    except Exception as e:
        log(f"Warning: failed to save clustermap: {e}")
