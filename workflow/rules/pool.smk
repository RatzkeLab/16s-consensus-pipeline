########################################
# 8) Build consensus database (FASTA)
########################################
rule build_consensus_db:
    input: db_inputs_from_checkpoint
    output:
        db = f"{OUT}/consensus_db.fasta"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.db})"
        cat {input} > {output.db}
        """