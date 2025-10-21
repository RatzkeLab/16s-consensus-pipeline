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
    Concatenate all multi-consensus sequences into a single database FASTA.
    Includes all cluster variants (e.g., sampleA_A, sampleA_B, etc.)
    """
    input:
        consensus_dirs=get_multi_consensus_dirs
    output:
        database=MULTI_DATABASE_FILE
    log:
        LOG_DIR / "pool" / "multi_db.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.database})"
        
        # Find all FASTA files in multi_consensus subdirectories
        find {MULTI_CONSENSUS_DIR} -name "*.fasta" -type f -exec cat {{}} + > {output.database}
        
        NUM_SEQS=$(grep -c "^>" {output.database} || echo "0")
        echo "Pooled $NUM_SEQS multi-consensus sequences into database" > {log}
        """