"""
Pool all consensus sequences into a single database FASTA file.
"""

rule pool_consensus:
    """
    Concatenate all consensus sequences into a single database FASTA.
    """
    input:
        consensus=get_consensus_files
    output:
        database=DATABASE_FILE
    log:
        LOG_DIR / "pool" / "consensus_db.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.database})"
        
        cat {input.consensus} > {output.database}
        
        NUM_SEQS=$(grep -c "^>" {output.database} || true)
        echo "Pooled $NUM_SEQS consensus sequences into database" > {log}
        """