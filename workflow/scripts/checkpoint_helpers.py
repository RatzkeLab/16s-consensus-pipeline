"""
Helper functions for handling Snakemake checkpoints and sample aggregation.

These functions depend on checkpoint outputs and are used to dynamically
determine which samples and files to process based on pipeline QC thresholds.
"""

from pathlib import Path


def read_passing_samples(checkpoint_output_path):
    """
    Load list of samples that passed a checkpoint.
    
    Args:
        checkpoint_output_path: Path to checkpoint passing samples file
    
    Returns:
        List of sample names that passed quality thresholds
    """
    passing = []
    with open(checkpoint_output_path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                passing.append(line)
    return passing


def get_passing_samples(checkpoints, checkpoint_output_path):
    """Load list of samples that passed the initial read count checkpoint."""
    return read_passing_samples(checkpoint_output_path)


def get_aligned_samples(checkpoints, checkpoint_output_path):
    """Load list of samples that passed the post-filter checkpoint."""
    return read_passing_samples(checkpoint_output_path)


def get_filtered_fastq_files(checkpoints, wildcards, filter_dir):
    """Get list of filtered FASTQ files for samples that passed checkpoint."""
    checkpoint_output = checkpoints.check_min_reads.get(**wildcards).output.passing
    passing = read_passing_samples(checkpoint_output)
    return [str(filter_dir / f"{sample}.fastq") for sample in passing]


def get_cluster_samples(checkpoints, wildcards, cluster_alignment_dir):
    """
    Get list of samples that have cluster alignments created.
    """
    checkpoint_output = checkpoints.check_min_reads_filtered.get(**wildcards).output.passing
    passing = read_passing_samples(checkpoint_output)
    
    cluster_samples = []
    for sample in passing:
        # Check if this sample has cluster alignments
        sample_cluster_dir = cluster_alignment_dir / sample
        if sample_cluster_dir.exists() and list(sample_cluster_dir.glob("*.fasta")):
            cluster_samples.append(sample)
    
    return cluster_samples


def get_naive_consensus_files(checkpoints, wildcards, naive_consensus_dir):
    """
    Get list of naive consensus FASTA files for samples that passed alignment.
    """
    checkpoint_output = checkpoints.check_min_reads_filtered.get(**wildcards).output.passing
    passing = read_passing_samples(checkpoint_output)
    return [str(naive_consensus_dir / f"{sample}.fasta") for sample in passing]


def get_cluster_fastqs_for_sample(checkpoints, wildcards):
    """
    Get list of cluster FASTQ files for a sample after split_reads checkpoint.
    
    Args:
        checkpoints: Snakemake checkpoints object
        wildcards: Snakemake wildcards object
    
    Returns:
        List of paths to cluster FASTQ files. Either [sample.fastq] if no 
        clusters detected, or [sample_A.fastq, sample_B.fastq, ...] if clusters detected
    """
    checkpoint_output = checkpoints.split_reads.get(**wildcards).output.outdir
    split_dir = Path(checkpoint_output)
    
    # Find all FASTQ files in the split directory
    fastq_files = list(split_dir.glob("*.fastq"))
    return [str(f) for f in sorted(fastq_files)]


def get_all_cluster_fastqs(checkpoints, wildcards):
    """
    Get all cluster FASTQ files for all passing samples.
    
    Args:
        checkpoints: Snakemake checkpoints object
        wildcards: Snakemake wildcards object
    
    Returns:
        List of paths to all cluster FASTQ files across all samples
    """
    checkpoint_output = checkpoints.check_min_reads_filtered.get(**wildcards).output.passing
    passing = read_passing_samples(checkpoint_output)
    all_fastqs = []
    
    for sample in passing:
        # Get cluster FASTQs for this sample
        sample_wildcards = type('obj', (object,), {'sample': sample})
        fastqs = get_cluster_fastqs_for_sample(checkpoints, sample_wildcards)
        all_fastqs.extend(fastqs)
    
    return all_fastqs


def get_all_cluster_consensus_files(checkpoints, wildcards, cluster_consensus_dir):
    """
    Get all cluster consensus files for all passing samples.
    
    Args:
        checkpoints: Snakemake checkpoints object
        wildcards: Snakemake wildcards object
        cluster_consensus_dir: Path to cluster consensus directory
    
    Returns:
        List of paths to all cluster consensus FASTA files
    """
    checkpoint_output = checkpoints.check_min_reads_filtered.get(**wildcards).output.passing
    passing = read_passing_samples(checkpoint_output)
    all_consensus = []
    
    for sample in passing:
        # Get split reads checkpoint output for this sample
        checkpoint_output = checkpoints.split_reads.get(sample=sample).output.outdir
        split_dir = Path(checkpoint_output)
        fastq_files = list(split_dir.glob("*.fastq"))
        
        # Generate consensus file paths for each cluster
        for fastq in sorted(fastq_files):
            cluster_name = fastq.stem
            consensus_file = cluster_consensus_dir / sample / f"{cluster_name}.fasta"
            all_consensus.append(str(consensus_file))
    
    return all_consensus
