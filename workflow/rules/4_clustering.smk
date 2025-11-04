"""
Cluster detection: detect and split multiple strains/variants within samples.

1. Generate per-read profiles from alignment based on variable positions
2. Detect subclusters using hierarchical clustering on profiles
3. Split reads by cluster assignment into separate FASTQ files (checkpoint)
4. Create new sequence alignments for each split cluster of reads

Pipeline position: FOURTH stage (parallel to naive_consensus.smk)
Upstream: 2_align.smk (rule initial_alignment)
Downstream: multi_consensus.smk (rule cluster_consensus)
"""

# ==================== Profile Generation ====================

rule generate_profiles:
    """
    Generate per-read profiles from alignment based on variable positions.
    
    Upstream: 2_align.smk (rule initial_alignment)
    Downstream: rule detect_clusters
    
    Identifies variable positions in the alignment and creates a profile
    for each read based on those positions. This step is separated from
    clustering to allow for:
    - Independent inspection of variable positions
    - Reuse of profiles for different clustering parameters
    - Easier debugging and QC
    
    Note: This rule always generates profiles, even if there are few
    variable positions. The decision whether to cluster is made by
    the downstream detect_clusters rule.
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        profile_tsv=PROFILE_DIR / "{sample}.tsv"
    params:
        min_agreement=MULTI_CONSENSUS_MIN_AGREEMENT,
        trim_bp=MULTI_CONSENSUS_TRIM_BP,
        auto_trim_flag=MULTI_CONSENSUS_AUTO_TRIM_FLAG,
        compress_gaps_flag=MULTI_CONSENSUS_COMPRESS_GAPS_FLAG,
    log:
        LOG_DIR / "generate_profiles" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/generate_profiles.py \
          {input.alignment} \
          {output.profile_tsv} \
          --min_agreement {params.min_agreement} \
          --trim_bp {params.trim_bp} \
          {params.auto_trim_flag} \
          {params.compress_gaps_flag} \
          2> {log}
        """


# ==================== Cluster Detection ====================

rule detect_clusters:
    """
    Detect if sample contains obvious subclusters from read profiles.
    
    Upstream: rule generate_profiles
    Downstream: checkpoint split_reads
    
    Performs hierarchical clustering on pre-generated read profiles.
    Makes all decisions about whether clustering should be performed:
    - Checks if there are enough variable positions
    - Checks if there are enough reads for clustering
    - Applies minimum cluster size thresholds
    - Limits maximum number of clusters
    
    Outputs cluster assignments if found, or marker file if single cluster.
    """
    input:
        profile_tsv=PROFILE_DIR / "{sample}.tsv"
    output:
        # Output directory will contain either:
        # - cluster_assignments.tsv (if clusters found)
        # - no_clusters.txt (if single cluster)
        outdir=directory(CLUSTER_DETECTION_DIR / "{sample}")
    params:
        min_cluster_size=MULTI_CONSENSUS_MIN_CLUSTER_SIZE,
        min_cluster_size_percent=MULTI_CONSENSUS_MIN_CLUSTER_SIZE_PERCENT,
        max_clusters=MULTI_CONSENSUS_MAX_CLUSTERS,
        min_variable_positions=MULTI_CONSENSUS_MIN_VARIABLE_POSITIONS,
    log:
        LOG_DIR / "detect_clusters" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/cluster_from_profiles.py \
          {input.profile_tsv} \
          {output.outdir} \
          --min_cluster_size {params.min_cluster_size} \
          --min_cluster_size_percent {params.min_cluster_size_percent} \
          --max_clusters {params.max_clusters} \
          --min_variable_positions {params.min_variable_positions} \
          2> {log}
        """


# ==================== Split Reads by Cluster ====================

checkpoint split_reads:
    """
    Split subsampled reads into cluster-specific FASTQ files.
    
    Upstream: rule subsample (in 2_alignment.smk) AND rule detect_clusters
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

