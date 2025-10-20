#!/usr/bin/env python3
import json
from pathlib import Path


def run():
    smk = globals().get("snakemake")
    if smk is None:
        raise SystemExit("This script is intended to be executed by Snakemake")

    filter_qc, sample_qc, align_qc, consensus_qc = smk.input
    out_html = smk.output["html"] if isinstance(smk.output, dict) else smk.output[0]

    def load(p):
        with open(p) as f:
            return json.load(f)

    fq = load(filter_qc)
    sq = load(sample_qc)
    aq = load(align_qc)
    cq = load(consensus_qc)

    Path(Path(out_html).parent).mkdir(parents=True, exist_ok=True)
    with open(out_html, 'w') as h:
        h.write("<html><head><title>Sample Report</title></head><body>\n")
        h.write("<h1>Per-sample QC report</h1>\n")
        h.write("<h2>Filtering</h2>\n")
        h.write(f"<pre>{json.dumps(fq, indent=2)}</pre>\n")
        h.write("<h2>Sampling</h2>\n")
        h.write(f"<pre>{json.dumps(sq, indent=2)}</pre>\n")
        h.write("<h2>Alignment</h2>\n")
        h.write(f"<pre>{json.dumps(aq, indent=2)}</pre>\n")
        h.write("<h2>Consensus</h2>\n")
        h.write(f"<pre>{json.dumps(cq, indent=2)}</pre>\n")
        h.write("</body></html>\n")


if __name__ == "__main__":
    run()
