"""
Generate per-cluster alignments for visual inspection.
"""

rule visualize_clusters:
    """
    Extract and save per-cluster alignments for quality checking.
    
    For each multi-consensus cluster, creates a separate alignment file
    containing only the reads assigned to that cluster. This allows visual
    inspection of cluster separation quality.
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta",
        cluster_dir=MULTI_CONSENSUS_DIR / "{sample}"
    output:
        # Output directory containing per-cluster alignment files
        outdir=directory(OUT_DIR / "cluster_alignments" / "{sample}")
    log:
        LOG_DIR / "visualize_clusters" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {output.outdir}
        
        python workflow/scripts/visualize_clusters.py \
          {input.alignment} \
          {input.cluster_dir} \
          {output.outdir} \
          {wildcards.sample} \
          2> {log}
        """
