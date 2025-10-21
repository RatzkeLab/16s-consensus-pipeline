
########################################
# 9) Simple HTML report (includes <50% sites summary)
########################################
rule report:
    input: report_inputs_from_checkpoint
    output:
        html = f"{OUT}/report/summary.html"
    log:
        "logs/report/build.log"
    script:
        "workflow/scripts/build_report.py"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.html})" "$(dirname {log})"
        python {script} > {output.html} 2> {log}
        """
