"""
Helper functions for common.smk configuration and setup.

This module contains utility functions used across the Snakemake pipeline
for handling sample names, building parameters, and determining input files.
"""

from pathlib import Path


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


def build_nanofilt_params(min_quality=0, min_length=0, max_length=0):
    """
    Build NanoFilt parameter string from parameters.
    
    Args:
        min_quality: Minimum average quality score (0 = skip)
        min_length: Minimum read length (0 = skip)
        max_length: Maximum read length (0 = skip)
        
    Returns:
        String of command-line arguments for NanoFilt
    """
    params = []
    
    if min_quality > 0:
        params.append(f"-q {min_quality}")
    
    if min_length > 0:
        params.append(f"-l {min_length}")
    
    if max_length > 0:
        params.append(f"--maxlength {max_length}")
    
    return " ".join(params)


def build_mafft_flags(algorithm, gap_open=0, gap_extend=0):
    """
    Build MAFFT flags from algorithm choice and gap penalties.
    
    Args:
        algorithm: Algorithm choice ("auto", "ginsi", or "default" or "")
        gap_open: Gap opening penalty (0 = use MAFFT defaults)
        gap_extend: Gap extension penalty (0 = use MAFFT defaults)
    
    Returns:
        String of MAFFT command-line flags
        
    Raises:
        ValueError: If algorithm is not recognized
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


def get_input_fastq_path(input_dir, sample):
    """
    Get the input FASTQ file path for a sample.
    
    Checks for both .fastq and .fastq.gz extensions.
    
    Args:
        input_dir: Path to input directory
        sample: Sample name
        
    Returns:
        Path to the input FASTQ file (as string)
    """
    input_path = Path(input_dir)
    
    # Check for .fastq first
    fastq_path = input_path / f"{sample}.fastq"
    if fastq_path.exists():
        return str(fastq_path)
    
    # Check for .fastq.gz
    fastq_gz_path = input_path / f"{sample}.fastq.gz"
    if fastq_gz_path.exists():
        return str(fastq_gz_path)
    
    # If neither exists, return .fastq (will trigger error in rule)
    return str(fastq_path)


def calculate_auto_trim(sequences):
    """
    Calculate auto-trim values based on longest leading/trailing gaps.
    
    Returns the position of the first non-gap character in the sequence with
    the longest leading gap, and the position of the last non-gap character
    in the sequence with the longest trailing gap.

    Used both in common_helpers.py and pairwise_distances.py.
    
    Args:
        sequences: Dictionary or iterable of sequences (aligned)
    
    Returns:
        tuple: (trim_start, trim_end) - number of positions to trim from start and end
    """
    if not sequences:
        return 0, 0
    
    # Handle both dict and list inputs
    seq_values = sequences.values() if isinstance(sequences, dict) else sequences
    
    max_leading_gaps = 0
    max_trailing_gaps = 0
    
    for seq in seq_values:
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

