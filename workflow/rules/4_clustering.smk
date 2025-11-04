"""
Cluster detection: detect and split multiple strains/variants within samples.

1. Detect subclusters in alignment (if multiple variants exist)
2. Split reads by cluster assignment into separate FASTQ files (checkpoint)
3. Realign each cluster separately for improved accuracy

Pipeline position: FOURTH stage (parallel to naive_consensus.smk)
Upstream: align.smk (rule align)
Downstream: multi_consensus.smk (rule cluster_consensus)
"""

# ==================== Cluster Detection ====================

rule detect_clusters:
    """
    Detect if sample contains obvious subclusters.
    
    Upstream: align.smk (rule align)
    Downstream: checkpoint split_reads
    
    Analyzes alignment to determine if multiple distinct variants exist.
    Outputs cluster assignments if found, or marker file if single cluster.
    The detection is based on:
    - Minimum agreement threshold for cluster separation
    - Minimum cluster size (absolute and percentage)
    - Maximum number of clusters to detect
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        # Output directory will contain either:
        # - cluster_assignments.tsv (if clusters found)
        # - no_clusters.txt (if single cluster)
        outdir=directory(CLUSTER_DETECTION_DIR / "{sample}")
    params:
        sample="{sample}",
        min_agreement=MULTI_CONSENSUS_MIN_AGREEMENT,
        min_cluster_size=MULTI_CONSENSUS_MIN_CLUSTER_SIZE,
        min_cluster_size_percent=MULTI_CONSENSUS_MIN_CLUSTER_SIZE_PERCENT,
        max_clusters=MULTI_CONSENSUS_MAX_CLUSTERS,
        min_variable_positions=MULTI_CONSENSUS_MIN_VARIABLE_POSITIONS,
        trim_bp=MULTI_CONSENSUS_TRIM_BP,
        auto_trim_flag=MULTI_CONSENSUS_AUTO_TRIM_FLAG,
        compress_gaps_flag=MULTI_CONSENSUS_COMPRESS_GAPS_FLAG,
    log:
        LOG_DIR / "detect_clusters" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/detect_clusters.py \
          {input.alignment} \
          {output.outdir} \
          {params.sample} \
          --min_agreement {params.min_agreement} \
          --min_cluster_size {params.min_cluster_size} \
          --min_cluster_size_percent {params.min_cluster_size_percent} \
          --max_clusters {params.max_clusters} \
          --min_variable_positions {params.min_variable_positions} \
          --trim_bp {params.trim_bp} \
          {params.auto_trim_flag} \
          {params.compress_gaps_flag} \
          2> {log}
        """


# ==================== Split Reads by Cluster ====================

checkpoint split_reads:
    """
    Split subsampled reads into cluster-specific FASTQ files.
    
    Upstream: align.smk (rule subsample) AND rule detect_clusters
    Downstream: rule realign_cluster (for each cluster detected)
    
    If clusters were detected, splits reads into separate files (sample_A.fastq, sample_B.fastq, etc.).
    If no clusters, copies the input file as-is (sample.fastq).
    
    This is a checkpoint because the number of output files varies dynamically
    depending on whether clusters were detected and how many.
    """
    input:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq",
        cluster_dir=CLUSTER_DETECTION_DIR / "{sample}"
    output:
        # Output directory will contain either:
        # - sample.fastq (if no clusters)
        # - sample_A.fastq, sample_B.fastq, etc. (if clusters found)
        outdir=directory(SPLIT_READS_DIR / "{sample}")
    params:
        sample="{sample}"
    log:
        LOG_DIR / "split_reads" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/split_reads.py \
          {input.fastq} \
          {input.cluster_dir} \
          {output.outdir} \
          {params.sample} \
          2> {log}
        """


# ==================== Realign Clusters ====================

rule realign_cluster:
    """
    Align reads for a specific cluster using MAFFT.
    
    Upstream: checkpoint split_reads (dynamically for each cluster)
    Downstream: multi_consensus.smk (rule cluster_consensus)
    
    After reads are split by cluster, each cluster is aligned separately
    to produce more accurate cluster-specific alignments. This improves
    consensus quality by avoiding cross-contamination between variants.
    
    Input: FASTQ file for a single cluster (or all reads if no clusters detected)
    Output: Aligned FASTA file for that cluster
    """
    input:
        fastq = SPLIT_READS_DIR / "{sample}" / "{cluster}.fastq"
    output:
        fasta = CLUSTER_ALIGNMENT_DIR / "{sample}" / "{cluster}.fasta"
    log:
        LOG_DIR / "realign_cluster" / "{sample}_{cluster}.log"
    conda:
        "../envs/align.yaml"
    params:
        mafft_flags=MAFFT_CLUSTER_ALIGN_FLAGS
    shell:
        """
        seqkit fq2fa {input.fastq} 2>> {log} | \
        mafft {params.mafft_flags} --thread 1 - > {output.fasta} 2>> {log}
        """

