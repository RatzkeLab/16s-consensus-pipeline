"""
Read filtering rules using NanoFilt.

This module filters ONT reads based on quality and length thresholds.
"""

# ==================== NanoFilt Filtering ====================

rule filter_reads:
    """
    Filter reads using NanoFilt based on quality and length thresholds.
    
    Applies the following filters (if configured > 0):
    - Minimum average Q-score
    - Minimum read length
    - Maximum read length
    
    Only processes samples that passed the initial read count check.
    """
    input:
        fastq=get_input_fastq
    output:
        fastq=FILTER_DIR / "{sample}.fastq"
    params:
        nanofilt_params=NANOFILT_PARAMS,
        sample="{sample}"
    log:
        LOG_DIR / "filter" / "{sample}.log"
    conda:
        "../envs/filter.yaml"
    threads: 1
    shell:
        """
        set -euo pipefail
        
        # Create output directory
        mkdir -p "$(dirname {output.fastq})" "$(dirname {log})"
        
        # Log filtering parameters
        echo "Filtering sample {params.sample}" > {log}
        echo "Parameters: {params.nanofilt_params}" >> {log}
        echo "" >> {log}
        
        # Apply NanoFilt (or pass through if no parameters)
        if [ -n "{params.nanofilt_params}" ]; then
            # Decompress if needed and filter
            if [[ {input.fastq} == *.gz ]]; then
                zcat {input.fastq} | NanoFilt {params.nanofilt_params} > {output.fastq} 2>> {log}
            else
                cat {input.fastq} | NanoFilt {params.nanofilt_params} > {output.fastq} 2>> {log}
            fi
            
            # Count filtered reads
            FILTERED_READS=$(awk 'END{{print NR/4}}' {output.fastq})
            echo "" >> {log}
            echo "Filtered reads: $FILTERED_READS" >> {log}
        else
            # No filtering - just copy/decompress
            if [[ {input.fastq} == *.gz ]]; then
                zcat {input.fastq} > {output.fastq}
            else
                cp {input.fastq} {output.fastq}
            fi
            echo "No filtering applied (pass-through)" >> {log}
        fi
        """