"""
Generate pipeline summary report.
"""

rule generate_summary:
    """
    Generate summary report showing sample attrition through pipeline stages.
    """
    input:
        read_summary = CHECK_DIR / "read_check_summary.tsv",
        filter_summary = CHECK_DIR / "filtered_check_summary.tsv",
        naive_db = NAIVE_DATABASE_FILE,
        multi_db = MULTI_DATABASE_FILE
    output:
        report = OUT_DIR / "pipeline_summary.md"
    conda:
        "../envs/qc.yaml"
    log:
        LOG_DIR / "summary" / "generate_summary.log"
    script:
        "../scripts/generate_summary.py"
