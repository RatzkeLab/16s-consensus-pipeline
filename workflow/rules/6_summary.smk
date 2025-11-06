"""
Summary report generation rules for the 16S consensus pipeline.

This module generates the final pipeline summary report, aggregating
statistics from all processing stages to show sample attrition and
final output counts.

Pipeline position: FINAL stage
Upstream: preprocessing.smk (checkpoints), naive_consensus.smk (rule pool_naive), 
         cluster_consensus.smk (rule pool_multi)
Downstream: None (terminal rule)
"""

rule pairwise_edit_distance:
    """
    Calculate pairwise edit distances between all consensus sequences.
    
    Upstream: cluster_consensus.smk (rule pool_multi)
    Downstream: None (analysis output)
    
    Computes edit distance (Levenshtein distance) between every pair of
    consensus sequences in the multi-consensus database to show similarity.
    """
    input:
        multi_db = MULTI_DATABASE_FILE
    output:
        distances = PAIRWISE_DISTANCE_FILE
    params:
        ignore_first_n_bp = PAIRWISE_DISTANCE_IGNORE_FIRST_N_BP,
        ignore_last_n_bp = PAIRWISE_DISTANCE_IGNORE_LAST_N_BP,
        auto_trim = PAIRWISE_DISTANCE_AUTO_TRIM
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "pairwise_distances.log"
    script:
        "../scripts/pairwise_distances.py"


rule pairwise_distance_heatmap:
    """
    Build a symmetric distance matrix from pairwise TSV, cluster it, and plot a heatmap.
    
    Upstream: pairwise_edit_distance
    Downstream: None (analysis output)
    """
    input:
        distances = PAIRWISE_DISTANCE_FILE
    output:
        matrix = PAIRWISE_DISTANCE_MATRIX_FILE,
        heatmap = PAIRWISE_DISTANCE_HEATMAP_FILE
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "pairwise_heatmap.log"
    script:
        "../scripts/distance_heatmap.py"


rule qc_alignment:
    """
    Create QC alignment of all consensus sequences using MAFFT.
    
    Upstream: cluster_consensus.smk (rule pool_multi)
    Downstream: rule profile_qc_alignment
    
    Aligns all consensus sequences together for user inspection and
    quick visual comparison of sequence similarity across samples.
    """
    input:
        multi_db = MULTI_DATABASE_FILE
    output:
        alignment = QC_ALIGNMENT_FILE
    conda:
        "../envs/align.yaml"
    log:
        LOG_DIR / "summary" / "qc_alignment.log"
    threads: 4
    params:
        mafft_flags=MAFFT_QC_ALIGN_FLAGS
    shell:
        """
        mafft {params.mafft_flags} --thread {threads} {input.multi_db} > {output.alignment} 2> {log}
        """


rule profile_qc_alignment:
    """
    Generate profiles from the final QC alignment of all consensus sequences.
    
    Upstream: rule qc_alignment
    Downstream: rule visualize_qc_alignment
    
    Creates profiles for each consensus sequence in the QC alignment to
    identify variable positions across all samples and prepare for
    clustering visualization.
    """
    input:
        alignment = QC_ALIGNMENT_FILE
    output:
        profile_tsv = QC_ALIGNMENT_PROFILES_DIR / "qc_alignment_profiles.tsv"
    params:
        min_minor_freq = QC_ALIGN_MIN_MINOR_FREQ,
        trim_bp = QC_ALIGN_TRIM_BP,
        auto_trim_flag = QC_ALIGN_AUTO_TRIM_FLAG,
        min_trim = QC_ALIGN_MIN_TRIM,
        max_trim = QC_ALIGN_MAX_TRIM,
        compress_gaps_flag = QC_ALIGN_COMPRESS_GAPS_FLAG,
    log:
        LOG_DIR / "summary" / "profile_qc_alignment.log"
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


rule visualize_qc_alignment:
    """
    Cluster and visualize all consensus sequences for final QC.
    
    Upstream: rule profile_qc_alignment
    Downstream: None (QC output)
    
    Generates hierarchical clustering and heatmap visualization of all
    consensus sequences to show:
    - Which samples are most similar
    - Whether distinct groups/species exist
    - Overall sequence diversity
    """
    input:
        profile_tsv = QC_ALIGNMENT_PROFILES_DIR / "qc_alignment_profiles.tsv"
    output:
        outdir = directory(QC_ALIGNMENT_CLUSTER_VIZ_DIR)
    params:
        min_cluster_size = QC_ALIGN_VIZ_MIN_CLUSTER_SIZE,  # For QC, allow single-sample clusters
        min_cluster_size_percent = QC_ALIGN_VIZ_MIN_CLUSTER_SIZE_PERCENT,
        max_clusters = QC_ALIGN_VIZ_MAX_CLUSTERS,  # Allow more clusters since this is cross-sample
        min_variable_positions = QC_ALIGN_VIZ_MIN_VARIABLE_POSITIONS,  # Lower threshold for final QC
    log:
        LOG_DIR / "summary" / "visualize_qc_alignment.log"
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


rule generate_summary:
    """
    Generate summary report showing sample attrition through pipeline stages.
    
    Upstream: All checkpoint summaries and pooled databases
    Downstream: None (terminal rule - produces final report)
    
    Aggregates statistics from:
    - Initial read count check (preprocessing)
    - Post-filter read count check (preprocessing)
    - Naive consensus database (naive_consensus)
    - Multi/cluster consensus database (cluster_consensus)
    
    Produces a markdown report showing how many samples passed each stage
    and how many consensus sequences were generated.
    """
    input:
        read_summary = CHECK_DIR / "read_check_summary.tsv",
        filter_summary = CHECK_DIR / "filtered_check_summary.tsv",
        naive_db = NAIVE_DATABASE_FILE,
        multi_db = MULTI_DATABASE_FILE
    output:
        report = PIPELINE_SUMMARY_FILE
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "generate_summary.log"
    script:
        "../scripts/generate_summary.py"
