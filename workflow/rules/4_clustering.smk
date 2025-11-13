"""
Cluster detection: detect and split multiple strains/variants within samples.

1. Generate per-read profiles from alignment based on variable positions
2. Detect subclusters using hierarchical clustering on profiles
3. Split reads by cluster assignment into separate FASTQ files (checkpoint)

Pipeline position: FOURTH stage (parallel to 3_naive_consensus.smk)
Upstream: 2_align.smk (rule initial_alignment)
Downstream: 5_cluster_consensus.smk (rule cluster_consensus)
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
        min_minor_freq=PROFILE_GEN_MIN_MINOR_FREQ,
        trim_bp=PROFILE_GEN_TRIM_BP,
        auto_trim=PROFILE_GEN_AUTO_TRIM,
        auto_trim_flag=PROFILE_GEN_AUTO_TRIM_FLAG,
        min_trim=PROFILE_GEN_MIN_TRIM,
        max_trim=PROFILE_GEN_MAX_TRIM,
        compress_gaps=PROFILE_GEN_COMPRESS_GAPS,
        compress_gaps_flag=PROFILE_GEN_COMPRESS_GAPS_FLAG,
    log:
        LOG_DIR / "generate_profiles" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/generate_profiles.py \
          {input.alignment} \
          {output.profile_tsv} \
          --min_minor_freq {params.min_minor_freq} \
          --trim_bp {params.trim_bp} \
          {params.auto_trim_flag} \
          --min_trim {params.min_trim} \
          --max_trim {params.max_trim} \
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
        # Visualization is generated in same directory with fixed name if enabled.
        outdir=directory(CLUSTER_DETECTION_DIR / "{sample}"),
        viz_figure=CLUSTER_DETECTION_DIR / "{sample}" / "profiles_dendrogram.png" if PROFILE_GEN_ENABLE_VIZ else []
    params:
        min_cluster_size=CLUSTERING_MIN_CLUSTER_SIZE,
        min_cluster_size_percent=CLUSTERING_MIN_CLUSTER_SIZE_PERCENT,
        max_clusters=CLUSTERING_MAX_CLUSTERS,
        min_variable_positions=CLUSTERING_MIN_VARIABLE_POSITIONS,
        viz_out_flag=CLUSTER_DETECTION_VIZ_FLAG,
    log:
        LOG_DIR / "detect_clusters" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
                python workflow/scripts/cluster_from_profiles.py \
                    {input.profile_tsv} \
                    {output.outdir} \
                    {params.viz_out_flag} \
                    --min_cluster_size {params.min_cluster_size} \
                    --min_cluster_size_percent {params.min_cluster_size_percent} \
                    --max_clusters {params.max_clusters} \
                    --min_variable_positions {params.min_variable_positions} \
                    2> {log}
        """

# ==================== Split Reads by Cluster ====================

checkpoint split_reads:
    """
    Split subsampled reads into cluster-specific FASTQ files (only if multiple clusters).
    
    Upstream: rule subsample (in 2_alignment.smk) AND rule detect_clusters
    Downstream: rule realign_cluster (for each cluster detected, only if multiple)
    
    If multiple clusters: splits reads into separate files (sample_A.fastq, sample_B.fastq, etc.)
    If single/no clusters: creates a marker file (no_split.txt) instead
    
    This is a checkpoint because the number of output files varies dynamically
    depending on whether multiple clusters were detected. When not split, downstream
    rules will use the naive consensus instead of recomputing.
    """
    input:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq",
        cluster_dir=CLUSTER_DETECTION_DIR / "{sample}"
    output:
        # Output directory will contain either:
        # - sample_A.fastq, sample_B.fastq, etc. (if multiple clusters)
        # - no_split.txt (if single cluster or no clusters)
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