#!/usr/bin/env python3
import json, csv
from pathlib import Path


def run():
    smk = globals().get("snakemake")
    if smk is None:
        raise SystemExit("This script is intended to be executed by Snakemake")

    inputs = smk.input
    out_tsv = smk.output["tsv"] if isinstance(smk.output, dict) else smk.output[0]

    rows = [(Path(p).stem, json.load(open(p))) for p in inputs]
    Path(Path(out_tsv).parent).mkdir(parents=True, exist_ok=True)
    with open(out_tsv, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(["sample", "length", "n_aligned", "low_conf_positions_count"]) 
        for stem, data in rows:
            w.writerow([stem, data.get("length"), data.get("n_aligned"), len(data.get("low_conf_positions", []))])


if __name__ == "__main__":
    run()
