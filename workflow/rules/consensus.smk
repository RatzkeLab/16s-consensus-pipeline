########################################
# 7) Plurality consensus; record <50% sites to CSV
########################################
rule naive_consensus:
    input:
        aln = f"{OUT}/aln/{{sample}}.aln.fasta"
    output:
        fa   = f"{OUT}/consensus/{{sample}}.fasta",
        csv  = f"{OUT}/consensus_stats/{{sample}}.csv"
    script:
        "workflow/scripts/naive_consensus.py"
    log:
        f"logs/consensus_call/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fa})" "$(dirname {output.csv})" "$(dirname {log})"
        python {script} > {output.fa} 2> {log}
        """