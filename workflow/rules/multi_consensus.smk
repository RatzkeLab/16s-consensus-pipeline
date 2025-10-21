"""
Multi-consensus generation for samples with multiple clusters.
"""

rule multi_consensus:
    """
    Generate multiple consensus sequences if sample contains multiple clusters.
    
    Analyzes alignment for positions with low agreement, creates per-read profiles,
    clusters reads based on Hamming distance, and generates consensus for each cluster.
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        # Output is a directory that will contain one or more FASTA files
        outdir=directory(MULTI_CONSENSUS_DIR / "{sample}")
    params:
        sample="{sample}",
        min_agreement=MULTI_CONSENSUS_MIN_AGREEMENT,
        min_cluster_size=MULTI_CONSENSUS_MIN_CLUSTER_SIZE,
        max_hamming=MULTI_CONSENSUS_MAX_HAMMING,
        min_prop=NAIVE_CONSENSUS_MIN_PROP
    log:
        LOG_DIR / "multi_consensus" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        
        python workflow/scripts/multi_consensus.py \
          {input.alignment} \
          {output.outdir} \
          {params.sample} \
          --min_agreement {params.min_agreement} \
          --min_cluster_size {params.min_cluster_size} \
          --max_hamming {params.max_hamming} \
          --min_consensus_prop {params.min_prop} \
          2> {log}
        """
