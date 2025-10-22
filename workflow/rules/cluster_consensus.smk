"""
Rule for generating consensus sequences for each cluster.

Each cluster gets its own naive consensus sequence, which will be
included in the multi_db.fasta output.
"""

rule cluster_consensus:
    """
    Generate naive consensus for a specific cluster.
    
    Input: Aligned FASTA for a cluster
    Output: Consensus FASTA for that cluster
    """
    input:
        fasta = CLUSTER_ALIGNMENT_DIR / "{sample}" / "{cluster}.fasta"
    output:
        consensus = CLUSTER_CONSENSUS_DIR / "{sample}" / "{cluster}.fasta"
    params:
        sample_name = lambda w: f"{w.sample}_{w.cluster}",
        variant_threshold = lambda w: config["variant_threshold"]
    log:
        LOG_DIR / "cluster_consensus" / "{sample}_{cluster}.log"
    conda:
        "../envs/align.yaml"
    script:
        "../scripts/naive_consensus.py"
