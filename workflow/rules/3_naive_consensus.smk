"""
Naive consensus: single consensus (1 sequence per sample without cluster detection).

This module handles the "naive" consensus path - generating a single consensus
sequence per sample without attempting to detect and separate multiple strains:
1. Generate consensus using winner-takes-all approach from alignment
2. Pool all naive consensuses into a single database FASTA

Pipeline position: THIRD stage (parallel to clustering.smk)
Upstream: 2_align.smk (rule initial_alignment)
Downstream: summary.smk (rule generate_summary)
"""

# ==================== Naive Consensus Generation ====================

rule naive_consensus:
    """
    Generate naive consensus sequence using winner-takes-all approach.
    
    Upstream: 2_align.smk (rule initial_alignment)
    Downstream: rule pool_naive
    
    Positions below threshold are marked with winner base and variants recorded.
    This assumes a single strain per sample and takes the most common base at
    each position (no cluster separation).
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        fasta=NAIVE_CONSENSUS_DIR / "{sample}.fasta",
        variants=NAIVE_CONSENSUS_DIR / "{sample}_variants.tsv"
    params:
        sample="{sample}",
        record_variants_below=NAIVE_CONSENSUS_RECORD_VARIANTS_BELOW
    log:
        LOG_DIR / "naive_consensus" / "{sample}.log"
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
          {params.record_variants_below} \
          2> {log}
        """


# ==================== Pool Naive Consensus Database ====================

rule pool_naive:
    """
    Concatenate all naive consensus sequences into a single database FASTA.
    
    Upstream: rule naive_consensus (all passing samples)
    Downstream: summary.smk (rule generate_summary)
    
    Aggregates all single-consensus sequences into one output file for
    downstream analysis (e.g., BLAST, taxonomy assignment).
    """
    input:
        consensus=get_naive_consensus_files
    output:
        database=NAIVE_DATABASE_FILE
    log:
        LOG_DIR / "pool" / "naive_db.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.database})"
        
        cat {input.consensus} > {output.database}
        
        NUM_SEQS=$(grep -c "^>" {output.database} || true)
        echo "Pooled $NUM_SEQS naive consensus sequences into database" > {log}
        """

