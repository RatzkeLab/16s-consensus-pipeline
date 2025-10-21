"""
Subsampling reads for alignment.
"""

rule subsample:
    """
    Randomly subsample reads to a maximum number for alignment.
    If sample has fewer reads than subsample_n, uses all reads.
    """
    input:
        fastq=FILTER_DIR / "{sample}.fastq"
    output:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq"
    params:
        n=SUBSAMPLE_N
    log:
        LOG_DIR / "subsample" / "{sample}.log"
    conda:
        "../envs/filter.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.fastq})"
        
        TOTAL_READS=$(awk 'END{{print NR/4}}' {input.fastq})
        
        if [ "$TOTAL_READS" -le "{params.n}" ] || [ "{params.n}" -eq 0 ]; then
            # Use all reads
            cp {input.fastq} {output.fastq}
            echo "Using all $TOTAL_READS reads (threshold: {params.n})" > {log}
        else
            # Randomly subsample
            seqtk sample -s 42 {input.fastq} {params.n} > {output.fastq} 2>> {log}
            echo "Subsampled {params.n} from $TOTAL_READS reads" >> {log}
        fi
        """
