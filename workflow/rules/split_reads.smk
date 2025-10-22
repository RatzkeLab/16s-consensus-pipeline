"""
Split subsampled reads by cluster assignment.
"""

checkpoint split_reads:
    """
    Split subsampled reads into cluster-specific FASTQ files.
    
    If clusters were detected, splits reads into separate files.
    If no clusters, copies the input file as-is.
    
    This is a checkpoint because the number of output files varies
    depending on whether clusters were detected.
    """
    input:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq",
        cluster_dir=OUT_DIR / "cluster_detection" / "{sample}"
    output:
        # Output directory will contain either:
        # - sample.fastq (if no clusters)
        # - sample_A.fastq, sample_B.fastq, etc. (if clusters found)
        outdir=directory(OUT_DIR / "split_reads" / "{sample}")
    params:
        sample="{sample}"
    log:
        LOG_DIR / "split_reads" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/split_reads.py \
          {input.fastq} \
          {input.cluster_dir} \
          {output.outdir} \
          {params.sample} \
          2> {log}
        """
