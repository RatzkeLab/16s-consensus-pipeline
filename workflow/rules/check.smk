"""
Quality control checkpoint rules for the 16S consensus pipeline.

This module contains checkpoint rules that filter samples based on quality metrics.
"""

# ==================== Per-Sample Read Count ====================

rule count_reads:
    """
    Count the number of reads in each sample's FASTQ file.
    
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
        
        # Write report
        echo -e "sample\treads" > {output.report}
        echo -e "{params.sample}\t$READ_COUNT" >> {output.report}
        
        # Log
        echo "Sample {params.sample}: $READ_COUNT reads" > {log}
        """


# ==================== Initial Read Count Checkpoint ====================

checkpoint check_min_reads:
    """
    Check if samples have minimum required number of reads.
    
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
        
        # Create output directory
        mkdir -p "$(dirname {output.passing})" "$(dirname {log})"
        
        # Initialize outputs
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
            
            # Add to summary
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
            echo "ERROR: No samples passed the minimum read threshold of {params.min_reads}" >> {log}
            exit 1
        fi
        """


# ==================== Post-Filter Read Count ====================

rule count_filtered_reads:
    """
    Count reads in filtered FASTQ files.
    """
    input:
        fastq=FILTER_DIR / "{sample}.fastq"
    output:
        report=CHECK_DIR / "{sample}_filtered_readcount.tsv"
    params:
        sample="{sample}"
    log:
        LOG_DIR / "check_filtered" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        READ_COUNT=$(awk 'END{{print NR/4}}' {input.fastq})
        echo -e "sample\treads" > {output.report}
        echo -e "{params.sample}\t$READ_COUNT" >> {output.report}
        echo "Sample {params.sample}: $READ_COUNT filtered reads" > {log}
        """


checkpoint check_min_reads_filtered:
    """
    Check if filtered samples have minimum required reads for alignment.
    """
    input:
        reports=get_filtered_fastq_files
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
        
        for fq in {input.reports}; do
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