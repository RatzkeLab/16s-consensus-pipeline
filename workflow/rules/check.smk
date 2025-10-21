########################################
# 1) Pre-filter read check (checkpoint)
########################################
checkpoint check_reads:
    input:
        fq = f"{DEMUX}/{{sample}}.fastq"
    output:
        tsv = f"{OUT}/checks/{{sample}}_pre.tsv"
    log:
        f"logs/check_reads/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
        READS=$(awk 'END{{print NR/4}}' {input.fq})
        printf "sample\treads\tthreshold\tstatus\n" > {output.tsv}
        STATUS="PASS"
        if [ "$READS" -lt "{MIN_READS_PRE}" ]; then STATUS="FAIL"; fi
        printf "%s\t%s\t%s\t%s\n" "{wildcards.sample}" "$READS" "{MIN_READS_PRE}" "$STATUS" >> {output.tsv}
        echo "Pre-filter check: $READS reads → $STATUS" > {log}
        """


########################################
# 3) Post-filter read check (checkpoint)
########################################
checkpoint postfilter_check:
    input:
        fq = f"{OUT}/filtered/{{sample}}.fastq"
    output:
        tsv = f"{OUT}/checks/{{sample}}_post.tsv"
    log:
        f"logs/postfilter_check/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.tsv})" "$(dirname {log})"
        READS=$(awk 'END{{print NR/4}}' {input.fq})
        printf "sample\treads\tthreshold\tstatus\n" > {output.tsv}
        STATUS="PASS"
        if [ "$READS" -lt "{MIN_READS_POST}" ]; then STATUS="FAIL"; fi
        printf "%s\t%s\t%s\t%s\n" "{wildcards.sample}" "$READS" "{MIN_READS_POST}" "$STATUS" >> {output.tsv}
        echo "Post-filter check: $READS reads → $STATUS" > {log}
        """