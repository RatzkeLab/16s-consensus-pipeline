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
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "pairwise_distances.log"
    script:
        "../scripts/pairwise_distances.py"


rule multi_align_consensus:
    """
    Create multi-alignment of all consensus sequences using MAFFT.
    
    Upstream: cluster_consensus.smk (rule pool_multi)
    Downstream: None (analysis output)
    
    Aligns all consensus sequences together for visualization and
    comparison of sequence variation across samples.
    """
    input:
        multi_db = MULTI_DATABASE_FILE
    output:
        alignment = MULTI_ALIGNMENT_FILE
    conda:
        "../envs/align.yaml"
    log:
        LOG_DIR / "summary" / "multi_align.log"
    threads: 4
    params:
        mafft_flags=MAFFT_MULTI_ALIGN_FLAGS
    shell:
        """
        mafft {params.mafft_flags} --thread {threads} {input.multi_db} > {output.alignment} 2> {log}
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
        report = OUT_DIR / "pipeline_summary.md"
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "generate_summary.log"
    script:
        "../scripts/generate_summary.py"
