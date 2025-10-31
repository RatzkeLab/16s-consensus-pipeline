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
    get_multi_consensus_dirs as _get_multi_consensus_dirs,
    get_cluster_fastqs_for_sample as _get_cluster_fastqs_for_sample,
    get_all_cluster_fastqs as _get_all_cluster_fastqs,
    get_all_cluster_consensus_files as _get_all_cluster_consensus_files,
)

# ==================== Configuration ====================

# Base directories
INPUT_DIR = Path(config.get("input_dir", "demo/merged"))
OUT_DIR = Path(config.get("outdir", "demo_output"))

# Structured output directories
CHECK_DIR = OUT_DIR / "checks"
FILTER_DIR = OUT_DIR / "filtered"
SUBSAMPLE_DIR = OUT_DIR / "subsampled"
ALIGNMENT_DIR = OUT_DIR / "alignment"
NAIVE_CONSENSUS_DIR = OUT_DIR / "naive_consensus"
CLUSTER_DETECTION_DIR = OUT_DIR / "cluster_detection"
SPLIT_READS_DIR = OUT_DIR / "split_reads"
CLUSTER_ALIGNMENT_DIR = OUT_DIR / "cluster_alignments"
CLUSTER_CONSENSUS_DIR = OUT_DIR / "cluster_consensus"
MULTI_CONSENSUS_DIR = OUT_DIR / "multi_consensus"
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

# ==================== Helper Functions ====================

def get_sample_names(input_dir):
    """
    Scan input directory for FASTQ files and return list of sample names.
    
    Supports both .fastq and .fastq.gz extensions.
    Sample name is the filename without extension(s).
    
    Args:
        input_dir: Path object or string pointing to input directory
        
    Returns:
        Sorted list of unique sample names
        
    Raises:
        ValueError: If no FASTQ files found in directory
    """
    input_path = Path(input_dir)
    sample_names = set()
    
    # Scan for .fastq files
    for fastq_file in input_path.glob("*.fastq"):
        sample_names.add(fastq_file.stem)
    
    # Scan for .fastq.gz files
    for fastq_gz in input_path.glob("*.fastq.gz"):
        # Remove .fastq.gz to get sample name
        sample_name = fastq_gz.name.replace(".fastq.gz", "")
        sample_names.add(sample_name)
    
    if not sample_names:
        raise ValueError(f"No FASTQ files found in {input_path}")
    
    return sorted(sample_names)


def build_nanofilt_params():
    """
    Build NanoFilt parameter string from config.
    
    Returns:
        String of command-line arguments for NanoFilt
    """
    params = []
    
    if NANOFILT_MIN_QUALITY > 0:
        params.append(f"-q {NANOFILT_MIN_QUALITY}")
    
    if NANOFILT_MIN_LENGTH > 0:
        params.append(f"-l {NANOFILT_MIN_LENGTH}")
    
    if NANOFILT_MAX_LENGTH > 0:
        params.append(f"--maxlength {NANOFILT_MAX_LENGTH}")
    
    return " ".join(params)

# MAFFT algorithm flags
def _build_mafft_flags(algorithm, gap_open=0, gap_extend=0):
    """
    Helper to build MAFFT flags from algorithm choice and gap penalties.
    
    Args:
        algorithm: Algorithm choice ("auto", "ginsi", or "default"/"")
        gap_open: Gap opening penalty (0 = use MAFFT defaults)
        gap_extend: Gap extension penalty (0 = use MAFFT defaults)
    
    Returns:
        String of MAFFT command-line flags
    """
    flags = []
    
    # Algorithm flags
    if algorithm == "default" or algorithm == "":
        pass  # No algorithm flag
    elif algorithm == "auto":
        flags.append("--auto")
    elif algorithm == "ginsi":
        flags.append("--globalpair --maxiterate 1000")
    else:
        raise ValueError(f"Unknown mafft_algorithm: {algorithm}. Must be 'auto' or 'ginsi'")
    
    # Gap penalty flags
    if gap_open > 0:
        flags.append(f"--op {gap_open}")
    if gap_extend > 0:
        flags.append(f"--ep {gap_extend}")
    
    return " ".join(flags)

def get_input_fastq(wildcards):
    """
    Get the input FASTQ file path for a sample.
    
    Checks for both .fastq and .fastq.gz extensions.
    
    Returns:
        Path to the input FASTQ file
    """
    sample = wildcards.sample
    
    # Check for .fastq first
    fastq_path = INPUT_DIR / f"{sample}.fastq"
    if fastq_path.exists():
        return str(fastq_path)
    
    # Check for .fastq.gz
    fastq_gz_path = INPUT_DIR / f"{sample}.fastq.gz"
    if fastq_gz_path.exists():
        return str(fastq_gz_path)
    
    # If neither exists, return .fastq (will trigger error in rule)
    return str(fastq_path)


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


def get_multi_consensus_dirs(wildcards):
    """Get multi-consensus output directories for all passing samples."""
    return _get_multi_consensus_dirs(checkpoints, wildcards, MULTI_CONSENSUS_DIR)


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

# Get all sample names from input directory
ALL_SAMPLES = get_sample_names(INPUT_DIR)

# Build NanoFilt parameters
NANOFILT_PARAMS = build_nanofilt_params()

# Build MAFFT flags for each alignment context
MAFFT_ALIGN_FLAGS = _build_mafft_flags(
    config.get("alignment", {}).get("mafft_algorithm", "default"),
    config.get("alignment", {}).get("gap_open", 0),
    config.get("alignment", {}).get("gap_extend", 0)
)
MAFFT_CLUSTER_ALIGN_FLAGS = _build_mafft_flags(
    config.get("cluster_alignment", {}).get("mafft_algorithm", "default"),
    config.get("cluster_alignment", {}).get("gap_open", 0),
    config.get("cluster_alignment", {}).get("gap_extend", 0)
)
MAFFT_MULTI_ALIGN_FLAGS = _build_mafft_flags(
    config.get("multi_alignment", {}).get("mafft_algorithm", "default"),
    config.get("multi_alignment", {}).get("gap_open", 0),
    config.get("multi_alignment", {}).get("gap_extend", 0)
)
