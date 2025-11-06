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
    get_all_cluster_alignment_viz_dirs as _get_all_cluster_alignment_viz_dirs,
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
CLUSTER_ALIGNMENT_PROFILES_DIR = OUT_DIR / "08b_cluster_alignment_profiles"
CLUSTER_ALIGNMENT_CLUSTER_VIZ_DIR = OUT_DIR / "08c_cluster_alignment_viz"
QC_DIR = OUT_DIR / "10_qc"
CLUSTER_DETECTION_VIZ_DIR = QC_DIR / "cluster_detection_viz"
QC_ALIGNMENT_PROFILES_DIR = QC_DIR / "qc_alignment_profiles"
QC_ALIGNMENT_CLUSTER_VIZ_DIR = QC_DIR / "qc_alignment_viz"
LOG_DIR = OUT_DIR / "logs"

# Thresholds and parameters
MIN_READS_INITIAL = config.get("min_reads_initial", 10)
MIN_READS_FILTERED = config.get("min_reads_filtered", 5)
SUBSAMPLE_N = config.get("subsample_n", 150)
SUBSAMPLE_SEED = config.get("subsample", {}).get("seed", 42)

# NanoFilt parameters
NANOFILT_MIN_QUALITY = config.get("filter", {}).get("min_avg_qscore", 10)
NANOFILT_MIN_LENGTH = config.get("filter", {}).get("min_length", 1300)
NANOFILT_MAX_LENGTH = config.get("filter", {}).get("max_length", 1700)
HEADCROP = config.get("filter", {}).get("headcrop", 0)
TAILCROP = config.get("filter", {}).get("tailcrop", 0)

# Parameters for Consensus Making
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

# Initial per-sample cluster-detection profile generation parameters
# (used by rule generate_profiles before clustering). Defaults fall back
# to multi_consensus values for convenience.
CLUSTER_DETECT_PROF_MIN_MINOR_FREQ = config.get("cluster_detection_profiles", {}).get(
    "min_minor_freq", 0.05
)
CLUSTER_DETECT_PROF_TRIM_BP = config.get("cluster_detection_profiles", {}).get(
    "trim_bp", 70
)
CLUSTER_DETECT_PROF_AUTO_TRIM = config.get("cluster_detection_profiles", {}).get(
    "auto_trim", True
)
CLUSTER_DETECT_PROF_AUTO_TRIM_FLAG = "--auto_trim" if CLUSTER_DETECT_PROF_AUTO_TRIM else ""
CLUSTER_DETECT_PROF_MIN_TRIM = config.get("cluster_detection_profiles", {}).get(
    "min_trim", 50
)
CLUSTER_DETECT_PROF_MAX_TRIM = config.get("cluster_detection_profiles", {}).get(
    "max_trim", 250
)
CLUSTER_DETECT_PROF_COMPRESS_GAPS = config.get("cluster_detection_profiles", {}).get(
    "compress_gaps", False
)
CLUSTER_DETECT_PROF_COMPRESS_GAPS_FLAG = "--compress_gaps" if CLUSTER_DETECT_PROF_COMPRESS_GAPS else ""
CLUSTER_DETECT_PROF_ENABLE_VIZ = config.get("cluster_detection_profiles", {}).get(
    "enable_viz", True
)

# Cluster-alignment QC profiling parameters (defaults fall back to cluster_detection_profiles)
CLUSTER_ALIGN_QC_MIN_MINOR_FREQ = config.get("cluster_alignment_qc", {}).get(
    "min_minor_freq", CLUSTER_DETECT_PROF_MIN_MINOR_FREQ
)
CLUSTER_ALIGN_QC_TRIM_BP = config.get("cluster_alignment_qc", {}).get(
    "trim_bp", CLUSTER_DETECT_PROF_TRIM_BP
)
CLUSTER_ALIGN_QC_AUTO_TRIM = config.get("cluster_alignment_qc", {}).get(
    "auto_trim", CLUSTER_DETECT_PROF_AUTO_TRIM
)
CLUSTER_ALIGN_QC_AUTO_TRIM_FLAG = "--auto_trim" if CLUSTER_ALIGN_QC_AUTO_TRIM else ""
CLUSTER_ALIGN_QC_MIN_TRIM = config.get("cluster_alignment_qc", {}).get(
    "min_trim", CLUSTER_DETECT_PROF_MIN_TRIM
)
CLUSTER_ALIGN_QC_MAX_TRIM = config.get("cluster_alignment_qc", {}).get(
    "max_trim", CLUSTER_DETECT_PROF_MAX_TRIM
)
CLUSTER_ALIGN_QC_COMPRESS_GAPS = config.get("cluster_alignment_qc", {}).get(
    "compress_gaps", CLUSTER_DETECT_PROF_COMPRESS_GAPS
)
CLUSTER_ALIGN_QC_COMPRESS_GAPS_FLAG = "--compress_gaps" if CLUSTER_ALIGN_QC_COMPRESS_GAPS else ""

# Cluster-alignment QC visualization clustering parameters
CLUSTER_ALIGN_QC_MIN_CLUSTER_SIZE = config.get("cluster_alignment_qc", {}).get(
    "min_cluster_size", MULTI_CONSENSUS_MIN_CLUSTER_SIZE
)
CLUSTER_ALIGN_QC_MIN_CLUSTER_SIZE_PERCENT = config.get("cluster_alignment_qc", {}).get(
    "min_cluster_size_percent", MULTI_CONSENSUS_MIN_CLUSTER_SIZE_PERCENT
)
CLUSTER_ALIGN_QC_MAX_CLUSTERS = config.get("cluster_alignment_qc", {}).get(
    "max_clusters", MULTI_CONSENSUS_MAX_CLUSTERS
)
CLUSTER_ALIGN_QC_MIN_VARIABLE_POSITIONS = config.get("cluster_alignment_qc", {}).get(
    "min_variable_positions", MULTI_CONSENSUS_MIN_VARIABLE_POSITIONS
)

# QC-alignment (all-consensus) profiling parameters (can be different defaults)
QC_ALIGN_MIN_MINOR_FREQ = config.get("qc_alignment_qc", {}).get(
    "min_minor_freq", CLUSTER_DETECT_PROF_MIN_MINOR_FREQ
)
QC_ALIGN_TRIM_BP = config.get("qc_alignment_qc", {}).get(
    "trim_bp", CLUSTER_DETECT_PROF_TRIM_BP
)
QC_ALIGN_AUTO_TRIM = config.get("qc_alignment_qc", {}).get(
    "auto_trim", CLUSTER_DETECT_PROF_AUTO_TRIM
)
QC_ALIGN_AUTO_TRIM_FLAG = "--auto_trim" if QC_ALIGN_AUTO_TRIM else ""
QC_ALIGN_MIN_TRIM = config.get("qc_alignment_qc", {}).get(
    "min_trim", CLUSTER_DETECT_PROF_MIN_TRIM
)
QC_ALIGN_MAX_TRIM = config.get("qc_alignment_qc", {}).get(
    "max_trim", CLUSTER_DETECT_PROF_MAX_TRIM
)
QC_ALIGN_COMPRESS_GAPS = config.get("qc_alignment_qc", {}).get(
    "compress_gaps", CLUSTER_DETECT_PROF_COMPRESS_GAPS
)
QC_ALIGN_COMPRESS_GAPS_FLAG = "--compress_gaps" if QC_ALIGN_COMPRESS_GAPS else ""

# QC-alignment visualization clustering parameters (more permissive defaults)
QC_ALIGN_VIZ_MIN_CLUSTER_SIZE = config.get("qc_alignment_qc", {}).get(
    "min_cluster_size", 1
)
QC_ALIGN_VIZ_MIN_CLUSTER_SIZE_PERCENT = config.get("qc_alignment_qc", {}).get(
    "min_cluster_size_percent", 0.0
)
QC_ALIGN_VIZ_MAX_CLUSTERS = config.get("qc_alignment_qc", {}).get(
    "max_clusters", 20
)
QC_ALIGN_VIZ_MIN_VARIABLE_POSITIONS = config.get("qc_alignment_qc", {}).get(
    "min_variable_positions", 1
)

# Pairwise distance parameters
PAIRWISE_DISTANCE_IGNORE_FIRST_N_BP = config.get("pairwise_distance", {}).get("ignore_first_n_bp", 70)
PAIRWISE_DISTANCE_IGNORE_LAST_N_BP = config.get("pairwise_distance", {}).get("ignore_last_n_bp", 70)
PAIRWISE_DISTANCE_AUTO_TRIM = config.get("pairwise_distance", {}).get("auto_trim", False)

# Final outputs
NAIVE_DATABASE_FILE = OUT_DIR / config.get("naive_database_filename", "naive_db.fasta")
MULTI_DATABASE_FILE = OUT_DIR / config.get("multi_database_filename", "multi_db.fasta")
PAIRWISE_DISTANCE_FILE = QC_DIR / config.get("pairwise_distance_filename", "pairwise_distances.tsv")
QC_ALIGNMENT_FILE = QC_DIR / config.get("qc_alignment_filename", "consensus_qc_alignment.fasta")
PAIRWISE_DISTANCE_MATRIX_FILE = QC_DIR / config.get("pairwise_distance_matrix_filename", "pairwise_distance_matrix.tsv")
PAIRWISE_DISTANCE_HEATMAP_FILE = QC_DIR / config.get("pairwise_distance_heatmap_filename", "pairwise_distance_heatmap.png")
PIPELINE_SUMMARY_FILE = QC_DIR / "pipeline_summary.md"

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
    """Get all cluster consensus files for all passing samples."""
    return _get_all_cluster_consensus_files(checkpoints, wildcards, CLUSTER_CONSENSUS_DIR)

def get_all_cluster_alignment_viz_dirs(wildcards):
    """Get all cluster alignment visualization directories for all passing samples."""
    return _get_all_cluster_alignment_viz_dirs(checkpoints, wildcards, CLUSTER_ALIGNMENT_CLUSTER_VIZ_DIR)

# ==================== Initialize ====================
# Dynamically generate command flags from config on initialization

# Get all sample names from input directory
ALL_SAMPLES = get_sample_names(INPUT_DIR)

NANOFILT_PARAMS = build_nanofilt_params(
    min_quality=NANOFILT_MIN_QUALITY,
    min_length=NANOFILT_MIN_LENGTH,
    max_length=NANOFILT_MAX_LENGTH
)

MAFFT_INITIAL_ALIGN_FLAGS = build_mafft_flags(
    config.get("initial_alignment", {}).get("mafft_algorithm", "default"),
    config.get("initial_alignment", {}).get("gap_open", 0),
    config.get("initial_alignment", {}).get("gap_extend", 0)
)
MAFFT_CLUSTER_ALIGN_FLAGS = build_mafft_flags(
    config.get("cluster_alignment", {}).get("mafft_algorithm", "default"),
    config.get("cluster_alignment", {}).get("gap_open", 0),
    config.get("cluster_alignment", {}).get("gap_extend", 0)
)
MAFFT_QC_ALIGN_FLAGS = build_mafft_flags(
    config.get("qc_alignment", {}).get("mafft_algorithm", "default"),
    config.get("qc_alignment", {}).get("gap_open", 0),
    config.get("qc_alignment", {}).get("gap_extend", 0)
)

