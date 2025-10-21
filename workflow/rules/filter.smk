########################################
# 2) Filtering / trimming (optional via config)
########################################
rule filter_reads:
    input:
        fq = f"{DEMUX}/{{sample}}.fastq"
    output:
        fq = f"{OUT}/filtered/{{sample}}.fastq"
    log:
        f"logs/filter_reads/{{sample}}.log"
    params:
        nanofilt_args = NANOFILT_ARGS
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.fq})" "$(dirname {log})"
        if [ -n "{params.nanofilt_args}" ]; then
          cat {input.fq} \
            | NanoFilt {params.nanofilt_args} \
            > {output.fq} 2> {log}
        else
          cat {input.fq} > {output.fq}
          echo "No filtering applied (pass-through)" > {log}
        fi
        """