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
CLUSTER_DETECTION_DIR = OUT_DIR / "05_cluster_detection"
SPLIT_READS_DIR = OUT_DIR / "06_split_reads"
CLUSTER_ALIGNMENT_DIR = OUT_DIR / "07_cluster_alignments"
CLUSTER_CONSENSUS_DIR = OUT_DIR / "08_cluster_consensus"
LOG_DIR = OUT_DIR / "logs"

# Thresholds and parameters
MIN_READS_INITIAL = config.get("min_reads_initial", 10)
MIN_READS_FILTERED = config.get("min_reads_filtered", 5)
SUBSAMPLE_N = config.get("subsample_n", 150)
# Deterministic subsampling seed (seqtk -s). Allows reproducible subsampling.
# Read from config.subsample.seed if provided, else defaults to 42.
SUBSAMPLE_SEED = config.get("subsample", {}).get("seed", 42)
NAIVE_CONSENSUS_MIN_PROP = config.get("naive_consensus", {}).get("min_consensus_proportion", 0.6)
# For multi-consensus, default to the naive value, or 0.6 if neither is set
MULTI_CONSENSUS_MIN_PROP = config.get("multi_consensus", {}).get(
    "min_consensus_proportion",
    NAIVE_CONSENSUS_MIN_PROP if 'NAIVE_CONSENSUS_MIN_PROP' in globals() else 0.6,
)
MULTI_CONSENSUS_MIN_AGREEMENT = config.get("multi_consensus", {}).get("min_agreement", 0.8)
MULTI_CONSENSUS_MIN_CLUSTER_SIZE = config.get("multi_consensus", {}).get("min_cluster_size", 5)
MULTI_CONSENSUS_MIN_CLUSTER_SIZE_PERCENT = config.get("multi_consensus", {}).get("min_cluster_size_percent", 0.0)
MULTI_CONSENSUS_MAX_CLUSTERS = config.get("multi_consensus", {}).get("max_clusters", 10)
MULTI_CONSENSUS_MIN_VARIABLE_POSITIONS = config.get("multi_consensus", {}).get("min_variable_positions", 3)
MULTI_CONSENSUS_TRIM_BP = config.get("multi_consensus", {}).get("trim_bp", 70)
MULTI_CONSENSUS_AUTO_TRIM = config.get("multi_consensus", {}).get("auto_trim", True)
MULTI_CONSENSUS_AUTO_TRIM_FLAG = "--auto_trim" if MULTI_CONSENSUS_AUTO_TRIM else ""
MULTI_CONSENSUS_COMPRESS_GAPS = config.get("multi_consensus", {}).get("compress_gaps", False)
MULTI_CONSENSUS_COMPRESS_GAPS_FLAG = "--compress_gaps" if MULTI_CONSENSUS_COMPRESS_GAPS else ""

# Final outputs
NAIVE_DATABASE_FILE = OUT_DIR / config.get("naive_database_filename", "naive_db.fasta")
MULTI_DATABASE_FILE = OUT_DIR / config.get("multi_database_filename", "multi_db.fasta")
PAIRWISE_DISTANCE_FILE = OUT_DIR / config.get("pairwise_distance_filename", "pairwise_distances.tsv")
MULTI_ALIGNMENT_FILE = OUT_DIR / config.get("multi_alignment_filename", "all_consensus_alignment.fasta")
PAIRWISE_DISTANCE_MATRIX_FILE = OUT_DIR / config.get("pairwise_distance_matrix_filename", "pairwise_distance_matrix.tsv")
PAIRWISE_DISTANCE_HEATMAP_FILE = OUT_DIR / config.get("pairwise_distance_heatmap_filename", "pairwise_distance_heatmap.png")

# Pairwise distance parameters
PAIRWISE_DISTANCE_IGNORE_FIRST_N_BP = config.get("pairwise_distance", {}).get("ignore_first_n_bp", 70)
PAIRWISE_DISTANCE_IGNORE_LAST_N_BP = config.get("pairwise_distance", {}).get("ignore_last_n_bp", 70)
PAIRWISE_DISTANCE_AUTO_TRIM = config.get("pairwise_distance", {}).get("auto_trim", False)

# NanoFilt parameters
NANOFILT_MIN_QUALITY = config.get("filter", {}).get("min_avg_qscore", 10)
NANOFILT_MIN_LENGTH = config.get("filter", {}).get("min_length", 1300)
NANOFILT_MAX_LENGTH = config.get("filter", {}).get("max_length", 1700)
HEADCROP = config.get("filter", {}).get("headcrop", 0)
TAILCROP = config.get("filter", {}).get("tailcrop", 0)

# ==================== Helper Functions (Wrappers) ====================
# See common_helpers.py for base functions
# See checkpoint_helpers.py for checkpoint-specific functions

def get_input_fastq(wildcards):
    return get_input_fastq_path(INPUT_DIR, wildcards.sample)

# Wrapper functions for checkpoint helpers
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
    """Get all cluster consensus files for all passing samples."""
    return _get_all_cluster_consensus_files(checkpoints, wildcards, CLUSTER_CONSENSUS_DIR)


# ==================== Initialize ====================
# Dynamically generate command flags from config on initialization

# Get all sample names from input directory
ALL_SAMPLES = get_sample_names(INPUT_DIR)

# Build NanoFilt parameters
NANOFILT_PARAMS = build_nanofilt_params(
    min_quality=NANOFILT_MIN_QUALITY,
    min_length=NANOFILT_MIN_LENGTH,
    max_length=NANOFILT_MAX_LENGTH
)

# Build MAFFT flags for each alignment context
MAFFT_ALIGN_FLAGS = build_mafft_flags(
    config.get("alignment", {}).get("mafft_algorithm", "default"),
    config.get("alignment", {}).get("gap_open", 0),
    config.get("alignment", {}).get("gap_extend", 0)
)
MAFFT_CLUSTER_ALIGN_FLAGS = build_mafft_flags(
    config.get("cluster_alignment", {}).get("mafft_algorithm", "default"),
    config.get("cluster_alignment", {}).get("gap_open", 0),
    config.get("cluster_alignment", {}).get("gap_extend", 0)
)
MAFFT_MULTI_ALIGN_FLAGS = build_mafft_flags(
    config.get("multi_alignment", {}).get("mafft_algorithm", "default"),
    config.get("multi_alignment", {}).get("gap_open", 0),
    config.get("multi_alignment", {}).get("gap_extend", 0)
)
