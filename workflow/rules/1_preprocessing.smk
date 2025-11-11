"""
Preprocessing: quality control and filtering.

This module handles the initial quality control and filtering stages of the pipeline:
1. Count raw reads per sample
2. Filter samples based on minimum read count threshold (checkpoint)
3. Apply quality and length filtering with NanoFilt
4. Re-count filtered reads and apply second threshold (checkpoint)

Pipeline position: FIRST stage
Upstream: Raw input FASTQ files
Downstream: subsample.smk (subsampling of filtered reads)
"""

# ==================== Per-Sample Read Count ====================

rule count_reads:
    """
    Count the number of reads in each sample's FASTQ file.
    
    Upstream: Raw input FASTQ files (configured in config.yaml)
    Downstream: checkpoint check_min_reads (aggregates all counts)
    
    Produces a TSV report with read count for the sample.
    """
    input:
        fastq=get_input_fastq
    output:
        report=CHECK_DIR / "{sample}_readcount.tsv"
    params:
        sample="{sample}"
    log:
        LOG_DIR / "check_reads" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        set -euo pipefail
        
        # Create output directories
        mkdir -p "$(dirname {output.report})" "$(dirname {log})"
        
        # Count reads (FASTQ has 4 lines per read)
        if [[ {input.fastq} == *.gz ]]; then
            READ_COUNT=$(zcat {input.fastq} | awk 'END{{print NR/4}}')
        else
            READ_COUNT=$(awk 'END{{print NR/4}}' {input.fastq})
        fi
        
        echo -e "sample\treads" > {output.report}
        echo -e "{params.sample}\t$READ_COUNT" >> {output.report}
        echo "Sample {params.sample}: $READ_COUNT reads" > {log}
        """


# ==================== Initial Read Count Checkpoint ====================

checkpoint check_min_reads:
    """
    Check if samples have minimum required number of reads.
    
    Upstream: rule count_reads (all samples)
    Downstream: rule filter_reads (only passing samples)
    
    This checkpoint evaluates all samples' read counts and creates:
    - A list of passing samples (for downstream processing)
    - A summary report of all samples and their status
    
    Samples that fail this check are excluded from further analysis.
    """
    input:
        reports=expand(
            str(CHECK_DIR / "{sample}_readcount.tsv"),
            sample=ALL_SAMPLES
        )
    output:
        passing=CHECK_DIR / "passing_samples.txt",
        summary=CHECK_DIR / "read_check_summary.tsv"
    params:
        min_reads=MIN_READS_INITIAL
    log:
        LOG_DIR / "check_reads" / "checkpoint.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.passing})" "$(dirname {log})"
        
        echo -e "sample\treads\tthreshold\tstatus" > {output.summary}
        > {output.passing}
        
        # Process each sample report
        for report in {input.reports}; do
            # Extract sample name and read count
            SAMPLE=$(tail -n 1 "$report" | cut -f1)
            READS=$(tail -n 1 "$report" | cut -f2)
            
            # Check threshold
            if [ "$READS" -ge "{params.min_reads}" ]; then
                STATUS="PASS"
                echo "$SAMPLE" >> {output.passing}
            else
                STATUS="FAIL"
            fi
            
            echo -e "$SAMPLE\t$READS\t{params.min_reads}\t$STATUS" >> {output.summary}
        done
        
        # Log summary statistics
        TOTAL=$(tail -n +2 {output.summary} | wc -l)
        PASSED=$(grep -c "PASS" {output.summary} || true)
        FAILED=$(grep -c "FAIL" {output.summary} || true)
        
        echo "Read count check summary:" > {log}
        echo "  Total samples: $TOTAL" >> {log}
        echo "  Passed: $PASSED" >> {log}
        echo "  Failed: $FAILED" >> {log}
        echo "" >> {log}
        echo "Passing samples written to: {output.passing}" >> {log}
        
        # Check if any samples passed
        if [ "$PASSED" -eq 0 ]; then
            echo "ERROR: No samples passed minimum read threshold of {params.min_reads}" >> {log}
            exit 1
        fi
        """


# ==================== NanoFilt Filtering ====================

rule filter_reads:
    """
    Filter reads using NanoFilt based on quality and length thresholds.
    
    Upstream: checkpoint check_min_reads (only processes passing samples)
    Downstream: checkpoint check_min_reads_filtered (counts filtered reads)
    
    Applies the following filters (if configured > 0):
    - Minimum average Q-score
    - Minimum read length
    - Maximum read length
    - Headcrop (trim from start)
    - Tailcrop (trim from end)
    
    Only processes samples that passed the initial read count check.
    """
    input:
        fastq=get_input_fastq
    output:
        fastq=temp(FILTER_DIR / "{sample}.fastq")
    params:
        nanofilt_params=NANOFILT_PARAMS,
        headcrop=FILTER_HEADCROP,
        tailcrop=FILTER_TAILCROP,
        sample="{sample}"
    log:
        LOG_DIR / "filter" / "{sample}.log"
    conda:
        "../envs/filter.yaml"
    threads: 1
    shell:
        """
        set -euo pipefail
        mkdir -p "$(dirname {output.fastq})" "$(dirname {log})"
        
        # Log filtering parameters
        echo "Filtering sample {params.sample}" > {log}
        echo "Parameters: {params.nanofilt_params}" >> {log}
        echo "Headcrop: {params.headcrop}" >> {log}
        echo "Tailcrop: {params.tailcrop}" >> {log}
        echo "" >> {log}
        
        # Apply NanoFilt (or pass through if no parameters)
        if [ -n "{params.nanofilt_params}" ]; then
            # Decompress if needed and filter
            if [[ {input.fastq} == *.gz ]]; then
                zcat {input.fastq} | NanoFilt {params.nanofilt_params} --headcrop {params.headcrop} --tailcrop {params.tailcrop} > {output.fastq} 2>> {log}
            else
                cat {input.fastq} | NanoFilt {params.nanofilt_params} --headcrop {params.headcrop} --tailcrop {params.tailcrop} > {output.fastq} 2>> {log}
            fi
            
            # Count filtered reads
            FILTERED_READS=$(awk 'END{{print NR/4}}' {output.fastq})
            echo "" >> {log}
            echo "Filtered reads: $FILTERED_READS" >> {log}
        else
            if [[ {input.fastq} == *.gz ]]; then
                zcat {input.fastq} > {output.fastq}
            else
                cp {input.fastq} {output.fastq}
            fi
            echo "No filtering applied (pass-through)" >> {log}
        fi
        """


# ==================== Post-Filter Checkpoint ====================

checkpoint check_min_reads_filtered:
    """
    Check if filtered samples have minimum required reads for alignment.
    
    Upstream: rule count_filtered_reads (all filtered samples)
    Downstream: subsample.smk (rule subsample - only passing filtered samples)
    
    This second checkpoint ensures samples still have enough reads after
    quality filtering to proceed with alignment and consensus calling.
    """
    input:
        fastqs=get_filtered_fastq_files
    output:
        passing=CHECK_DIR / "passing_filtered_samples.txt",
        summary=CHECK_DIR / "filtered_check_summary.tsv"
    params:
        min_reads=MIN_READS_FILTERED
    log:
        LOG_DIR / "check_filtered" / "checkpoint.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.passing})"
        echo -e "sample\treads\tthreshold\tstatus" > {output.summary}
        > {output.passing}
        
        for fq in {input.fastqs}; do
            SAMPLE=$(basename "$fq" .fastq)
            READS=$(awk 'END{{print NR/4}}' "$fq")
            
            if [ "$READS" -ge "{params.min_reads}" ]; then
                STATUS="PASS"
                echo "$SAMPLE" >> {output.passing}
            else
                STATUS="FAIL"
            fi
            
            echo -e "$SAMPLE\t$READS\t{params.min_reads}\t$STATUS" >> {output.summary}
        done
        
        PASSED=$(grep -c "PASS" {output.summary} || true)
        echo "Filtered samples passed: $PASSED" > {log}
        
        if [ "$PASSED" -eq 0 ]; then
            echo "ERROR: No samples passed post-filter threshold of {params.min_reads}" >> {log}
            exit 1
        fi
        """

