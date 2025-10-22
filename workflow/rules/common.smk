from pathlib import Path

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
MULTI_CONSENSUS_DIR = OUT_DIR / "multi_consensus"
LOG_DIR = OUT_DIR / "logs"

# Thresholds and parameters
MIN_READS_INITIAL = config.get("min_reads_initial", 10)
MIN_READS_FILTERED = config.get("min_reads_filtered", 5)
SUBSAMPLE_N = config.get("subsample_n", 150)
NAIVE_CONSENSUS_MIN_PROP = config.get("naive_consensus", {}).get("min_consensus_proportion", 0.6)
MULTI_CONSENSUS_MIN_AGREEMENT = config.get("multi_consensus", {}).get("min_agreement", 0.8)
MULTI_CONSENSUS_MIN_CLUSTER_SIZE = config.get("multi_consensus", {}).get("min_cluster_size", 5)
MULTI_CONSENSUS_MAX_CLUSTERS = config.get("multi_consensus", {}).get("max_clusters", 10)

# Final outputs
NAIVE_DATABASE_FILE = OUT_DIR / config.get("naive_database_filename", "naive_db.fasta")
MULTI_DATABASE_FILE = OUT_DIR / config.get("multi_database_filename", "multi_db.fasta")

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


def get_passing_samples(wildcards):
    """
    Load list of samples that passed the initial read count checkpoint.
    
    Reads from the checkpoint output file containing passing sample names.
    
    Returns:
        List of sample names that passed quality thresholds
    """
    checkpoint_output = checkpoints.check_min_reads.get(**wildcards).output.passing
    
    passing = []
    with open(checkpoint_output) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                passing.append(line)
    
    return passing


def get_filtered_fastq_files(wildcards):
    """
    Get list of filtered FASTQ files for samples that passed checkpoint.
    
    This function is used in aggregate rules that need all filtered files.
    
    Returns:
        List of paths to filtered FASTQ files
    """
    passing = get_passing_samples(wildcards)
    return [str(FILTER_DIR / f"{sample}.fastq") for sample in passing]


def get_aligned_samples(wildcards):
    """
    Load list of samples that passed the post-filter checkpoint.
    
    Returns:
        List of sample names that passed filtered read count threshold
    """
    checkpoint_output = checkpoints.check_min_reads_filtered.get(**wildcards).output.passing
    
    passing = []
    with open(checkpoint_output) as f:
        for line in f:
            line = line.strip()
            if line:
                passing.append(line)
    
    return passing


def get_naive_consensus_files(wildcards):
    """
    Get list of naive consensus FASTA files for samples that passed alignment.
    """
    passing = get_aligned_samples(wildcards)
    return [str(NAIVE_CONSENSUS_DIR / f"{sample}.fasta") for sample in passing]


def get_multi_consensus_dirs(wildcards):
    """
    Get list of multi-consensus output directories for all passing samples.
    Each sample has its own subdirectory that may contain multiple cluster files.
    """
    passing = get_aligned_samples(wildcards)
    return [str(MULTI_CONSENSUS_DIR / sample) for sample in passing]


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


# ==================== Initialize ====================

# Get all sample names from input directory
ALL_SAMPLES = get_sample_names(INPUT_DIR)

# Build NanoFilt parameters
NANOFILT_PARAMS = build_nanofilt_params()
