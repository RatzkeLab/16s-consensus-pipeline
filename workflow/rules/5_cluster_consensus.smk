"""
Cluster-based consensus generation rules for the 16S consensus pipeline.

This module handles consensus calling for samples that have been split into
multiple clusters (representing different strains/variants):
1. Re-align the split reads for each cluster
2. Generate consensus for each individual cluster
3. Pool all cluster consensuses into a single database FASTA

Pipeline position: FIFTH stage (parallel to naive_consensus.smk)
Upstream: 4_clustering.smk (rule realign_cluster)
Downstream: 6_summary.smk (rule generate_summary)
"""

# ==================== Realign Clusters ====================

def get_cluster_fastq(wildcards):
    """
    Dynamic input function for realign_cluster.
    Only collects FASTQ files from split_reads directory (ignores no_split samples).
    """
    # Trigger the checkpoint
    checkpoints.split_reads.get(sample=wildcards.sample)
    
    # Return the FASTQ file path for this specific cluster
    return str(SPLIT_READS_DIR / wildcards.sample / f"{wildcards.cluster}.fastq")


rule realign_cluster:
    """
    Align reads for a specific cluster using MAFFT.
    
    Upstream: checkpoint split_reads (only for samples with multiple clusters)
    Downstream: rule cluster_consensus
    
    After reads are split by cluster, each cluster is aligned separately
    to produce more accurate cluster-specific alignments. This improves
    consensus quality by avoiding cross-contamination between variants.
    
    Note: This rule only runs for samples that were actually split into 
    multiple clusters. Single-cluster samples use naive consensus instead.
    
    Input: FASTQ file for a single cluster
    Output: Aligned FASTA file for that cluster
    """
    input:
        fastq = get_cluster_fastq
    output:
        fasta = CLUSTER_ALIGNMENT_DIR / "{sample}" / "{cluster}.fasta"
    log:
        LOG_DIR / "realign_cluster" / "{sample}_{cluster}.log"
    conda:
        "../envs/align.yaml"
    params:
        mafft_flags=MAFFT_SECONDARY_ALIGN_FLAGS
    shell:
        """
        seqkit fq2fa {input.fastq} 2>> {log} | \
        mafft {params.mafft_flags} --thread 1 - > {output.fasta} 2>> {log}
        """



# ==================== Per-Cluster Consensus Generation ====================

rule cluster_consensus:
    """
    Generate naive consensus for a specific cluster.
    
    Upstream: clustering.smk (rule realign_cluster)
    Downstream: rule pool_multi
    
    Takes the cluster-specific alignment and generates a consensus sequence
    using the same winner-takes-all approach as naive consensus, but applied
    to each cluster independently for improved accuracy.
    
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
        record_variants_below = CLUSTER_CONSENSUS_RECORD_VARIANTS_BELOW
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
          {params.record_variants_below} \
          2> {log}
        """


# ==================== Pool Cluster Consensus Database ====================

def get_multi_consensus_inputs(wildcards):
    """
    Get all consensus files for multi database.
    - For split samples: cluster consensus files
    - For non-split samples: naive consensus files
    """
    # Get split samples (have multiple clusters)
    split_samples = get_split_samples(wildcards)
    split_consensus = []
    for sample in split_samples:
        # Trigger checkpoint
        checkpoints.split_reads.get(sample=sample)
        split_dir = SPLIT_READS_DIR / sample
        fastq_files = list(split_dir.glob("*.fastq"))
        for fastq in sorted(fastq_files):
            cluster_name = fastq.stem
            consensus_file = CLUSTER_CONSENSUS_DIR / sample / f"{cluster_name}.fasta"
            split_consensus.append(str(consensus_file))
    
    # Get non-split samples (single cluster - use naive consensus)
    no_split_samples = get_no_split_samples(wildcards)
    naive_consensus = [str(NAIVE_CONSENSUS_DIR / f"{sample}.fasta") for sample in no_split_samples]
    
    return split_consensus + naive_consensus


rule pool_multi:
    """
    Concatenate all consensus sequences into a single database FASTA.
    
    Upstream: rule cluster_consensus (for split samples), rule naive_consensus (for non-split samples)
    Downstream: summary.smk (rule generate_summary)
    
    For each sample:
    - If multiple clusters detected: includes all cluster consensuses (e.g., sample_A, sample_B)
    - If single/no clusters: includes the naive consensus for that sample (avoids redundant computation)
    
    This uses the split_reads checkpoint to dynamically determine which samples
    were split vs. which should use naive consensus.
    """
    input:
        consensus_files = get_multi_consensus_inputs
    output:
        database = MULTI_DATABASE_FILE
    log:
        LOG_DIR / "pool" / "multi_db.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.database})"
        
        # Concatenate all consensus files
        cat {input.consensus_files} > {output.database}
        
        NUM_SEQS=$(grep -c "^>" {output.database} || echo "0")
        echo "Pooled $NUM_SEQS consensus sequences into database" > {log}
        """


# (Per-cluster alignment QC via per-read profiling and clustering removed as requested.)

