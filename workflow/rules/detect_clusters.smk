"""
Detect subclusters in samples for multi-consensus generation.
"""

rule detect_clusters:
    """
    Detect if sample contains obvious subclusters.
    
    Analyzes alignment to determine if multiple distinct variants exist.
    Outputs cluster assignments if found, or marker file if single cluster.
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        # Output directory will contain either:
        # - cluster_assignments.tsv (if clusters found)
        # - no_clusters.txt (if single cluster)
        outdir=directory(OUT_DIR / "cluster_detection" / "{sample}")
    params:
        sample="{sample}",
        min_agreement=MULTI_CONSENSUS_MIN_AGREEMENT,
        min_cluster_size=MULTI_CONSENSUS_MIN_CLUSTER_SIZE,
        min_cluster_size_percent=MULTI_CONSENSUS_MIN_CLUSTER_SIZE_PERCENT,
        max_clusters=MULTI_CONSENSUS_MAX_CLUSTERS
    log:
        LOG_DIR / "detect_clusters" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/detect_clusters.py \
          {input.alignment} \
          {output.outdir} \
          {params.sample} \
          --min_agreement {params.min_agreement} \
          --min_cluster_size {params.min_cluster_size} \
          --min_cluster_size_percent {params.min_cluster_size_percent} \
          --max_clusters {params.max_clusters} \
          2> {log}
        """
