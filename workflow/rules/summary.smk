"""
Generate pipeline summary report.
"""

rule generate_summary:
    """
    Create a summary report showing sample attrition through the pipeline.
    """
    input:
        initial_check=CHECK_DIR / "read_check_summary.tsv",
        filtered_check=CHECK_DIR / "filtered_check_summary.tsv",
        database=DATABASE_FILE
    output:
        summary=OUT_DIR / "pipeline_summary.md"
    log:
        LOG_DIR / "summary" / "generate.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        python workflow/scripts/generate_summary.py \
          {input.initial_check} \
          {input.filtered_check} \
          {output.summary} \
          2> {log}
        
        echo "Summary report generated at {output.summary}" >> {log}
        """
