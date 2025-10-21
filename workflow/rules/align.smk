########################################
# 6) Multiple sequence alignment (MAFFT)
########################################
rule align_reads:
    input:
        fq = f"{OUT}/subsampled/{{sample}}.fastq"
    output:
        aln = f"{OUT}/aln/{{sample}}.aln.fasta"
    threads: MAFFT_THREADS
    log:
        f"logs/align_reads/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.aln})" "$(dirname {log})"
        seqtk seq -A {input.fq} \
          | mafft --thread {threads} --auto - \
          > {output.aln} 2> {log}
        """
