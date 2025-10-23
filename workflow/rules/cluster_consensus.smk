"""
Naive consensus sequence generation for each cluster alignment.

This mirrors the naive_consensus rule but runs per-cluster alignment output.
"""

rule cluster_consensus:
    """
    Generate naive consensus for a specific cluster.
    
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
