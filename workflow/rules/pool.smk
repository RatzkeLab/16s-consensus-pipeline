"""
Pool consensus sequences into database FASTA files.
"""

rule pool_naive:
    """
    Concatenate all naive consensus sequences into a single database FASTA.
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


rule pool_multi:
    """
    Concatenate all cluster consensus sequences into a single database FASTA.
    
    For each sample:
    - If clusters detected: includes all cluster consensuses (e.g., sampleA_A, sampleA_B)
    - If no clusters: includes the naive consensus for that sample
    
    This uses the split_reads checkpoint to dynamically discover cluster files.
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
