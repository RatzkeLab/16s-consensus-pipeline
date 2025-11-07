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

rule realign_cluster:
    """
    Align reads for a specific cluster using MAFFT.
    
    Upstream: checkpoint split_reads (dynamically for each cluster)
    Downstream: multi_consensus.smk (rule cluster_consensus)
    
    After reads are split by cluster, each cluster is aligned separately
    to produce more accurate cluster-specific alignments. This improves
    consensus quality by avoiding cross-contamination between variants.
    
    Input: FASTQ file for a single cluster (or all reads if no clusters detected)
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
    params:
        mafft_flags=MAFFT_CLUSTER_ALIGN_FLAGS
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


# ==================== Pool Cluster Consensus Database ====================

rule pool_multi:
    """
    Concatenate all cluster consensus sequences into a single database FASTA.
    
    Upstream: rule cluster_consensus (all clusters from all samples)
    Downstream: summary.smk (rule generate_summary)
    
    For each sample:
    - If clusters detected: includes all cluster consensuses (e.g., sampleA_A, sampleA_B)
    - If no clusters: includes the naive consensus for that sample
    
    This uses the split_reads checkpoint to dynamically discover cluster files.
    Aggregates all split consensus sequences into one output file for
    downstream analysis (e.g., BLAST, taxonomy assignment).
    """
    input:
        # Use checkpoint to trigger split_reads for all samples
        split_dirs = lambda w: [checkpoints.split_reads.get(sample=s).output.outdir 
                                for s in get_aligned_samples(w)],
        # Cluster consensus files for all clusters
        cluster_consensus = lambda w: get_all_cluster_consensus_files(w)
    output:
        database = MULTI_DATABASE_FILE
    log:
        LOG_DIR / "pool" / "multi_db.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.database})"
        
        # Concatenate all cluster consensus files
        cat {input.cluster_consensus} > {output.database}
        
        NUM_SEQS=$(grep -c "^>" {output.database} || echo "0")
        echo "Pooled $NUM_SEQS cluster consensus sequences into database" > {log}
        """


# (Per-cluster alignment QC via per-read profiling and clustering removed as requested.)

