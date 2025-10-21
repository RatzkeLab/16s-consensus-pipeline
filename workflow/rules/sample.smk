########################################
# 5) Subsample reads (~N) — only for PASSED samples
########################################
rule sample_reads:
    input:
        fq = f"{OUT}/filtered/{{sample}}.fastq"
    output:
        fq = f"{OUT}/subsampled/{{sample}}.fastq"
    log:
        f"logs/sample_reads/{{sample}}.log"
    params:
        n = SUBS_N,
        seed = SUBS_SEED
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fq})" "$(dirname {log})"
        seqtk sample -s {params.seed} {input.fq} {params.n} > {output.fq} 2> {log}
        """