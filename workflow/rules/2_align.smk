"""
Initial alignment: subsampling and initial alignment of filtered reads.

1. Subsample reads to a maximum number (optional, for efficiency)
2. Perform the initial multiple sequence alignment using MAFFT

This initial alignment output is used by both the naive consensus path and the
cluster detection/multi-consensus path.

Pipeline position: SECOND stage
Upstream: preprocessing.smk (checkpoint check_min_reads_filtered)
Downstream: naive_consensus.smk (rule naive_consensus) AND clustering.smk (rule detect_clusters)
"""

# ==================== Subsampling ====================

rule subsample:
    """
    Randomly subsample reads to a maximum number for alignment.
    
    Upstream: preprocessing.smk (checkpoint check_min_reads_filtered)
    Downstream: rule initial_alignment
    
    If sample has fewer reads than subsample_n, uses all reads.
    Subsampling is deterministic (uses configured seed for reproducibility).
    This step is optional but recommended for samples with very high read
    counts to reduce computational cost while maintaining representative coverage.
    """
    input:
        fastq=FILTER_DIR / "{sample}.fastq"
    output:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq"
    params:
        n=SUBSAMPLE_N,
        seed=SUBSAMPLE_SEED
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
            # Randomly subsample (deterministic with configured seed)
            seqtk sample -s {params.seed} {input.fastq} {params.n} > {output.fastq} 2>> {log}
            echo "Subsampled {params.n} from $TOTAL_READS reads" >> {log}
        fi
        """


# ==================== Multiple Sequence Alignment ====================

rule initial_alignment:
    """
    Create multiple sequence alignment per sample using MAFFT.
    
    Upstream: rule subsample
    Downstream: naive_consensus.smk (rule naive_consensus) AND clustering.smk (rule detect_clusters)
    
    This alignment is used for both pipeline paths:
    - Naive consensus: direct consensus calling from alignment
    - Multi-consensus: cluster detection followed by per-cluster realignment and consensus
    """
    input:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq"
    output:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    log:
        LOG_DIR / "initial_alignment" / "{sample}.log"
    conda:
        "../envs/align.yaml"
    threads: 4
    params:
        initial_mafft_flags=MAFFT_INITIAL_ALIGN_FLAGS
    shell:
        """
        mkdir -p "$(dirname {output.alignment})"
        
        # Convert FASTQ to FASTA and align with MAFFT
        awk 'NR%4==1{{print ">"substr($0,2)}} NR%4==2{{print}}' {input.fastq} \
          | mafft --thread {threads} {params.initial_mafft_flags} - \
          > {output.alignment} 2> {log}
        """

