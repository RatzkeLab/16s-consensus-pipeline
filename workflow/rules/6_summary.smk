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
