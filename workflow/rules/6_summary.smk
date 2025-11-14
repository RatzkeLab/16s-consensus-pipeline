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

rule profile_from_qc_alignment:
    """
    Generate a profile (variant table) from the QC alignment with a generous MAF threshold (1%).
    Includes nearly all columns/positions in the alignment.
    """
    input:
        alignment = QC_ALIGNMENT_FILE
    output:
        profile = QC_DIR / "qc_alignment_profile.tsv"
    params:
        maf = 0.01
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "qc_alignment_profile.log"
    shell:
        """
        python workflow/scripts/generate_profiles.py \
          {input.alignment} \
          {output.profile} \
          --min_minor_freq {params.maf} \
          --trim_bp 0 \
          --compress_gaps \
          2> {log}
        """

rule cluster_from_qc_profile:
    """
    Cluster the QC alignment profile and generate a clustering heatmap.
    Uses a very inclusive profile (MAF >1%) to cluster all consensus sequences.
    """
    input:
        profile = QC_DIR / "qc_alignment_profile.tsv"
    output:
        heatmap = QC_DIR / "qc_profile_clustering" / "qc_profile_clustering_heatmap.png"
    params:
        outdir = QC_DIR / "qc_profile_clustering"
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "qc_profile_clustering.log"
    shell:
        """
        mkdir -p {params.outdir}
    python workflow/scripts/cluster_from_profiles.py {input.profile} {params.outdir} --viz_out $(basename {output.heatmap}) 2> {log}
        """


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


rule global_consensus:
    """
    Generate a global consensus sequence from all consensus sequences.
    
    Upstream: cluster_consensus.smk (rule pool_multi)
    Downstream: None (analysis output)
    
    Creates a single consensus sequence representing the most common
    base at each position across all consensus sequences in the multi
    database. Useful for overall QC and reference.
    """
    input:
        global_alignment = QC_ALIGNMENT_FILE
    output:
        fasta = QC_DIR / "global_consensus.fasta",
        variants = QC_DIR / "global_consensus_variants.tsv"
    params:
        sample = "global_consensus",
        record_variants_below = GLOBAL_CONSENSUS_RECORD_VARIANTS_BELOW
    log:
        LOG_DIR / "summary" / "global_consensus.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.fasta})"
        
        python workflow/scripts/naive_consensus.py \
          {input.global_alignment} \
          {output.fasta} \
          {output.variants} \
          {params.sample} \
          {params.record_variants_below} \
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
