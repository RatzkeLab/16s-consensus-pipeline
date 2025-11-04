"""
Cluster-based consensus generation rules for the 16S consensus pipeline.

This module handles consensus calling for samples that have been split into
multiple clusters (representing different strains/variants):
1. Generate consensus for each individual cluster
2. Pool all cluster consensuses into a single database FASTA

Pipeline position: FIFTH stage (parallel to naive_consensus.smk)
Upstream: clustering.smk (rule realign_cluster)
Downstream: summary.smk (rule generate_summary)
"""

# ==================== Per-Cluster Consensus Generation ====================

rule cluster_consensus:
    """
    Generate naive consensus for a specific cluster.
    
    Upstream: clustering.smk (rule realign_cluster)
    Downstream: rule pool_multi
    
    Takes the cluster-specific alignment and generates a consensus sequence
    using the same winner-takes-all approach as naive consensus, but applied
    to each cluster independently for improved accuracy.
    
    Input: Aligned FASTA for a cluster
    Output: Consensus FASTA (+ variants TSV) for that cluster
    """
    input:
        alignment = CLUSTER_ALIGNMENT_DIR / "{sample}" / "{cluster}.fasta"
    output:
        fasta = CLUSTER_CONSENSUS_DIR / "{sample}" / "{cluster}.fasta",
        variants = CLUSTER_CONSENSUS_DIR / "{sample}" / "{cluster}_variants.tsv"
    params:
        sample = lambda w: f"{w.sample}_{w.cluster}",
        min_prop = MULTI_CONSENSUS_MIN_PROP
    log:
        LOG_DIR / "cluster_consensus" / "{sample}_{cluster}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.fasta})"
        
        python workflow/scripts/naive_consensus.py \
          {input.alignment} \
          {output.fasta} \
          {output.variants} \
          {params.sample} \
          {params.min_prop} \
          2> {log}
        """


# ==================== Pool Cluster Consensus Database ====================

rule pool_multi:
    """
    Concatenate all cluster consensus sequences into a single database FASTA.
    
    Upstream: rule cluster_consensus (all clusters from all samples)
    Downstream: summary.smk (rule generate_summary)
    
    For each sample:
    - If clusters detected: includes all cluster consensuses (e.g., sampleA_A, sampleA_B)
    - If no clusters: includes the naive consensus for that sample
    
    This uses the split_reads checkpoint to dynamically discover cluster files.
    Aggregates all split consensus sequences into one output file for
    downstream analysis (e.g., BLAST, taxonomy assignment).
    """
    input:
        # Use checkpoint to trigger split_reads for all samples
        split_dirs = lambda w: [checkpoints.split_reads.get(sample=s).output.outdir 
                                for s in get_aligned_samples(w)],
        # Cluster consensus files for all clusters
        cluster_consensus = lambda w: get_all_cluster_consensus_files(w)
    output:
        database = MULTI_DATABASE_FILE
    log:
        LOG_DIR / "pool" / "multi_db.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.database})"
        
        # Concatenate all cluster consensus files
        cat {input.cluster_consensus} > {output.database}
        
        NUM_SEQS=$(grep -c "^>" {output.database} || echo "0")
        echo "Pooled $NUM_SEQS cluster consensus sequences into database" > {log}
        """


# ==================== Cluster Alignment QC (Profile & Visualize) ====================

rule profile_cluster_alignment:
    """
    Generate read profiles from cluster alignments for QC visualization.
    
    Upstream: clustering.smk (rule realign_cluster)
    Downstream: rule visualize_cluster_alignment
    
    Creates profiles for each read in the cluster-specific alignment to
    validate that the cluster alignment is internally consistent and that
    reads within the cluster are indeed similar.
    """
    input:
        alignment = CLUSTER_ALIGNMENT_DIR / "{sample}" / "{cluster}.fasta"
    output:
        profile_tsv = CLUSTER_ALIGNMENT_PROFILES_DIR / "{sample}" / "{cluster}.tsv"
    params:
        min_agreement = CLUSTER_ALIGN_QC_MIN_AGREEMENT,
        trim_bp = CLUSTER_ALIGN_QC_TRIM_BP,
        auto_trim_flag = CLUSTER_ALIGN_QC_AUTO_TRIM_FLAG,
        compress_gaps_flag = CLUSTER_ALIGN_QC_COMPRESS_GAPS_FLAG,
    log:
        LOG_DIR / "profile_cluster_alignment" / "{sample}_{cluster}.log"
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


rule visualize_cluster_alignment:
    """
    Cluster and visualize reads within each cluster alignment for QC.
    
    Upstream: rule profile_cluster_alignment
    Downstream: None (QC output)
    
    Generates clustering visualization to validate that:
    - Reads within a cluster are homogeneous
    - No sub-clusters were missed
    - Alignment quality is good
    """
    input:
        profile_tsv = CLUSTER_ALIGNMENT_PROFILES_DIR / "{sample}" / "{cluster}.tsv"
    output:
        outdir = directory(CLUSTER_ALIGNMENT_CLUSTER_VIZ_DIR / "{sample}" / "{cluster}")
    params:
        min_cluster_size = CLUSTER_ALIGN_QC_MIN_CLUSTER_SIZE,
        min_cluster_size_percent = CLUSTER_ALIGN_QC_MIN_CLUSTER_SIZE_PERCENT,
        max_clusters = CLUSTER_ALIGN_QC_MAX_CLUSTERS,
        min_variable_positions = CLUSTER_ALIGN_QC_MIN_VARIABLE_POSITIONS,
    log:
        LOG_DIR / "visualize_cluster_alignment" / "{sample}_{cluster}.log"
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


# ==================== Aggregate Cluster Alignment QC ====================

rule aggregate_cluster_alignment_qc:
    """
    Aggregate rule to trigger all cluster alignment visualizations.
    
    This rule ensures that visualizations are generated for all cluster
    alignments across all samples. It doesn't produce output itself but
    depends on all the visualization directories.
    """
    input:
        lambda w: get_all_cluster_alignment_viz_dirs(w)
    output:
        touch(QC_DIR / "cluster_alignment_qc_complete.txt")
    log:
        LOG_DIR / "aggregate" / "cluster_alignment_qc.log"
    shell:
        """
        echo "Generated QC visualizations for $(echo {input} | wc -w) cluster alignments" > {log}
        """

