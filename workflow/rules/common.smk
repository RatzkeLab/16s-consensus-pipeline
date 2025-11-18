from pathlib import Path
import sys

# Add scripts directory to path for imports
sys.path.insert(0, str(Path(workflow.basedir) / "scripts"))

from checkpoint_helpers import (
    get_passing_samples as _get_passing_samples,
    get_aligned_samples as _get_aligned_samples,
    get_filtered_fastq_files as _get_filtered_fastq_files,
    get_cluster_samples as _get_cluster_samples,
    get_naive_consensus_files as _get_naive_consensus_files,
    get_cluster_fastqs_for_sample as _get_cluster_fastqs_for_sample,
    get_all_cluster_fastqs as _get_all_cluster_fastqs,
    get_all_cluster_consensus_files as _get_all_cluster_consensus_files,
    get_split_samples as _get_split_samples,
    get_no_split_samples as _get_no_split_samples,
)

from common_helpers import (
    get_sample_names,
    build_nanofilt_params,
    build_mafft_flags,
    get_input_fastq_path,
)

# ==================== Configuration ====================

# Base directories
INPUT_DIR = Path(config.get("input_dir", "demo/merged"))
OUT_DIR = Path(config.get("outdir", "demo_output"))

# Structured output directories
CHECK_DIR = OUT_DIR / "00_checks"
FILTER_DIR = OUT_DIR / "01_filtered"
SUBSAMPLE_DIR = OUT_DIR / "02_subsampled"
ALIGNMENT_DIR = OUT_DIR / "03_alignment"
NAIVE_CONSENSUS_DIR = OUT_DIR / "04_naive_consensus"
PROFILE_DIR = OUT_DIR / "05_read_profiles"
CLUSTER_DETECTION_DIR = OUT_DIR / "06_cluster_detection"
SPLIT_READS_DIR = OUT_DIR / "07_split_reads"
CLUSTER_ALIGNMENT_DIR = OUT_DIR / "08_cluster_alignments"
CLUSTER_CONSENSUS_DIR = OUT_DIR / "09_cluster_consensus"
QC_DIR = OUT_DIR / "10_qc"
LOG_DIR = OUT_DIR / "logs"

# Thresholds and parameters
MIN_READS_INITIAL = config.get("min_reads_initial", 10)
MIN_READS_FILTERED = config.get("min_reads_filtered", 5)
SUBSAMPLE_N = config.get("subsample_n", 150)
SUBSAMPLE_SEED = config.get("subsample", {}).get("seed", 42)

# Filter parameters
FILTER_HEADCROP = config.get("filter", {}).get("headcrop", 0)
FILTER_TAILCROP = config.get("filter", {}).get("tailcrop", 0)
FILTER_MIN_AVG_QSCORE = config.get("filter", {}).get("min_avg_qscore", 10)
FILTER_MIN_LENGTH = config.get("filter", {}).get("min_length", 1300)
FILTER_MAX_LENGTH = config.get("filter", {}).get("max_length", 1700)

# Naive consensus parameters
NAIVE_CONSENSUS_RECORD_VARIANTS_BELOW = config.get("naive_consensus", {}).get("record_variants_below", 0.6)

# Initial alignment parameters
INITIAL_ALIGNMENT_MAFFT_ALGORITHM = config.get("initial_alignment", {}).get("mafft_algorithm", "ginsi")
INITIAL_ALIGNMENT_GAP_OPEN = config.get("initial_alignment", {}).get("gap_open", 0)
INITIAL_ALIGNMENT_GAP_EXTEND = config.get("initial_alignment", {}).get("gap_extend", 0)

# Per-read profile generation parameters (for cluster detection)
PROFILE_GEN_MIN_MINOR_FREQ = config.get("per_read_profile_generation", {}).get("min_minor_freq", 0.05)
PROFILE_GEN_TRIM_BP = config.get("per_read_profile_generation", {}).get("trim_bp", 70)
PROFILE_GEN_AUTO_TRIM = config.get("per_read_profile_generation", {}).get("auto_trim", True)
PROFILE_GEN_MIN_TRIM = config.get("per_read_profile_generation", {}).get("min_trim", 50)
PROFILE_GEN_MAX_TRIM = config.get("per_read_profile_generation", {}).get("max_trim", 250)
PROFILE_GEN_COMPRESS_GAPS = config.get("per_read_profile_generation", {}).get("compress_gaps", False)
PROFILE_GEN_ENABLE_VIZ = config.get("per_read_profile_generation", {}).get("enable_viz", True)

# Clustering algorithm parameters (per-sample clustering)
CLUSTERING_MIN_VARIABLE_POSITIONS = config.get("clustering_algorithm", {}).get("min_variable_positions", 3)
CLUSTERING_MIN_CLUSTER_SIZE = config.get("clustering_algorithm", {}).get("min_cluster_size", 5)
CLUSTERING_MIN_CLUSTER_SIZE_PERCENT = config.get("clustering_algorithm", {}).get("min_cluster_size_percent", 5.0)
CLUSTERING_MAX_CLUSTERS = config.get("clustering_algorithm", {}).get("max_clusters", 10)
CLUSTERING_METHOD = config.get("clustering_algorithm", {}).get("clustering_method", "hdbscan")

# QC clustering parameters (final consensus sequence clustering)
QC_CLUSTERING_METHOD = config.get("qc_clustering", {}).get("clustering_method", "hdbscan")
QC_CLUSTERING_MAX_CLUSTERS = config.get("qc_clustering", {}).get("max_clusters", 10)

# Secondary alignment (cluster realignment) parameters
SECONDARY_ALIGNMENT_MAFFT_ALGORITHM = config.get("secondary_alignment", {}).get("mafft_algorithm", "ginsi")
SECONDARY_ALIGNMENT_GAP_OPEN = config.get("secondary_alignment", {}).get("gap_open", 0)
SECONDARY_ALIGNMENT_GAP_EXTEND = config.get("secondary_alignment", {}).get("gap_extend", 0)

# Cluster consensus parameters
CLUSTER_CONSENSUS_RECORD_VARIANTS_BELOW = config.get("cluster_consensus", {}).get("record_variants_below", 0.6)

# Pairwise distance parameters
PAIRWISE_DISTANCE_INCLUDE_NAIVE = config.get("pairwise_distance", {}).get("include_naive_consensus", True)
PAIRWISE_DISTANCE_AUTO_TRIM = config.get("pairwise_distance", {}).get("auto_trim", True)
PAIRWISE_DISTANCE_IGNORE_FIRST_N_BP = config.get("pairwise_distance", {}).get("ignore_first_n_bp", 70)
PAIRWISE_DISTANCE_IGNORE_LAST_N_BP = config.get("pairwise_distance", {}).get("ignore_last_n_bp", 70)

# QC alignment parameters
QC_ALIGNMENT_MAFFT_ALGORITHM = config.get("qc_alignment", {}).get("mafft_algorithm", "default")
QC_ALIGNMENT_GAP_OPEN = config.get("qc_alignment", {}).get("gap_open", 0)
QC_ALIGNMENT_GAP_EXTEND = config.get("qc_alignment", {}).get("gap_extend", 0)

# Global consensus parameters
GLOBAL_CONSENSUS_RECORD_VARIANTS_BELOW = config.get("global_consensus", {}).get("record_variants_below", 0.4)

# Final outputs
NAIVE_DATABASE_FILE = OUT_DIR / config.get("naive_database_filename", "naive_db.fasta")
MULTI_DATABASE_FILE = OUT_DIR / config.get("multi_database_filename", "multi_db.fasta")
PAIRWISE_DISTANCE_FILE = QC_DIR / "pairwise_distances.tsv"
QC_ALIGNMENT_FILE = QC_DIR / "consensus_qc_alignment.fasta"
PAIRWISE_DISTANCE_MATRIX_FILE = QC_DIR / "pairwise_distance_matrix.tsv"
PIPELINE_SUMMARY_FILE = QC_DIR / "pipeline_summary.md"
QC_PROFILE_CLUSTERING_HEATMAP = QC_DIR / "qc_profile_clustering" / "qc_profile_clustering_heatmap.png"
QC_PROFILE_DISTANCE_HEATMAP = QC_DIR / "qc_profile_clustering" / "distance_heatmap.png"

# ==================== Helper Functions (Wrappers) ====================
# See common_helpers.py for base functions
# See checkpoint_helpers.py for checkpoint-specific functions

def get_input_fastq(wildcards):
    return get_input_fastq_path(INPUT_DIR, wildcards.sample)

def get_passing_samples(wildcards):
    """Get samples that passed initial read count checkpoint."""
    checkpoint_output = checkpoints.check_min_reads.get(**wildcards).output.passing
    return _get_passing_samples(checkpoints, checkpoint_output)

def get_filtered_fastq_files(wildcards):
    """Get filtered FASTQ files for samples that passed checkpoint."""
    return _get_filtered_fastq_files(checkpoints, wildcards, FILTER_DIR)

def get_aligned_samples(wildcards):
    """Get samples that passed post-filter checkpoint."""
    checkpoint_output = checkpoints.check_min_reads_filtered.get(**wildcards).output.passing
    return _get_aligned_samples(checkpoints, checkpoint_output)

def get_cluster_samples(wildcards):
    """Get samples that have cluster alignments."""
    return _get_cluster_samples(checkpoints, wildcards, CLUSTER_ALIGNMENT_DIR)

def get_naive_consensus_files(wildcards):
    """Get naive consensus FASTA files for aligned samples."""
    return _get_naive_consensus_files(checkpoints, wildcards, NAIVE_CONSENSUS_DIR)

def get_cluster_fastqs_for_sample(wildcards):
    """Get cluster FASTQ files for a sample after split_reads checkpoint."""
    return _get_cluster_fastqs_for_sample(checkpoints, wildcards)

def get_all_cluster_fastqs(wildcards):
    """Get all cluster FASTQ files for all passing samples."""
    return _get_all_cluster_fastqs(checkpoints, wildcards)

def get_all_cluster_consensus_files(wildcards):
    """Get all cluster consensus files for all passing samples that were split."""
    return _get_all_cluster_consensus_files(checkpoints, wildcards, CLUSTER_CONSENSUS_DIR)

def get_split_samples(wildcards):
    """Get samples that were split into multiple clusters."""
    return _get_split_samples(checkpoints, wildcards, SPLIT_READS_DIR)

def get_no_split_samples(wildcards):
    """Get samples that were NOT split (single/no cluster)."""
    return _get_no_split_samples(checkpoints, wildcards, SPLIT_READS_DIR)

# ==================== Initialize ====================
# Dynamically generate command flags from config on initialization

# Get all sample names from input directory
ALL_SAMPLES = get_sample_names(INPUT_DIR)

# Build NanoFilt parameters string
NANOFILT_PARAMS = build_nanofilt_params(
    min_quality=FILTER_MIN_AVG_QSCORE,
    min_length=FILTER_MIN_LENGTH,
    max_length=FILTER_MAX_LENGTH
)

# Build MAFFT flags for different alignment stages
MAFFT_INITIAL_ALIGN_FLAGS = build_mafft_flags(
    INITIAL_ALIGNMENT_MAFFT_ALGORITHM,
    INITIAL_ALIGNMENT_GAP_OPEN,
    INITIAL_ALIGNMENT_GAP_EXTEND
)

MAFFT_SECONDARY_ALIGN_FLAGS = build_mafft_flags(
    SECONDARY_ALIGNMENT_MAFFT_ALGORITHM,
    SECONDARY_ALIGNMENT_GAP_OPEN,
    SECONDARY_ALIGNMENT_GAP_EXTEND
)

MAFFT_QC_ALIGN_FLAGS = build_mafft_flags(
    QC_ALIGNMENT_MAFFT_ALGORITHM,
    QC_ALIGNMENT_GAP_OPEN,
    QC_ALIGNMENT_GAP_EXTEND
)

# ==================== Flag Builders ====================
# Build CLI flags for boolean options so rules can inject them directly in shell blocks

def _bool_flag(enabled, flag: str):
    """Return flag string if enabled else empty string."""
    return flag if enabled else ""

def _viz_flag(enabled, output_path):
    """Return --viz_out flag with path if enabled, else empty string."""
    return f"--viz_out {output_path}" if enabled else ""

# Per-read profile generation flags
PROFILE_GEN_AUTO_TRIM_FLAG = _bool_flag(PROFILE_GEN_AUTO_TRIM, "--auto_trim")
PROFILE_GEN_COMPRESS_GAPS_FLAG = _bool_flag(PROFILE_GEN_COMPRESS_GAPS, "--compress_gaps")

# Cluster detection viz flag - built once since it's the same for all samples
def _get_cluster_viz_flag(filename: str):
    """Return a --viz_out flag using only a filename (no directory logic) if viz enabled.

    Caller is responsible for ensuring the working directory / output directory.
    This keeps rule params free of inline logic.
    """
    return f"--viz_out {filename}" if PROFILE_GEN_ENABLE_VIZ else ""

# Pre-built cluster detection visualization flag (fixed filename expected by downstream logic)
CLUSTER_DETECTION_VIZ_FLAG = _get_cluster_viz_flag("profiles_dendrogram.png")

