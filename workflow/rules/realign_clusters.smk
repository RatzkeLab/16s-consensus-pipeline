"""
Rule for realigning cluster-specific reads.

After reads are split by cluster, each cluster is aligned separately
to produce more accurate cluster-specific alignments.
"""

rule realign_cluster:
    """
    Align reads for a specific cluster using MAFFT.
    
    Input: FASTQ file for a single cluster (or all reads if no clusters)
    Output: Aligned FASTA file for that cluster
    """
    input:
        fastq = SPLIT_READS_DIR / "{sample}" / "{cluster}.fastq"
    output:
        fasta = CLUSTER_ALIGNMENT_DIR / "{sample}" / "{cluster}.fasta"
    log:
        LOG_DIR / "realign_cluster" / "{sample}_{cluster}.log"
    conda:
        "../envs/align.yaml"
    shell:
        """
        seqkit fq2fa {input.fastq} 2>> {log} | \
        mafft --auto --thread 1 - > {output.fasta} 2>> {log}
        """
