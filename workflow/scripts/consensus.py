#!/usr/bin/env python3
import json
from pathlib import Path


def read_fasta(path):
    seqs = []
    with open(path) as f:
        name = None
        parts = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    seqs.append((name, ''.join(parts)))
                name = line[1:].strip()
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seqs.append((name, ''.join(parts)))
    return seqs


def consensus_wta_with_conflict(aligned, threshold=0.5):
    # winner-takes-all per column, also record positions where winner <= threshold
    if not aligned:
        return "", []
    seqs = [s for _, s in aligned]
    length = max(len(s) for s in seqs)
    # pad with gaps
    padded = [s.ljust(length, '-') for s in seqs]
    cols = list(zip(*padded))
    cons = []
    low_conf_positions = []
    for i, col in enumerate(cols, start=1):
        counts = {}
        total_non_gap = 0
        for c in col:
            if c == '-':
                continue
            c = c.upper()
            counts[c] = counts.get(c, 0) + 1
            total_non_gap += 1
        if total_non_gap == 0:
            cons.append('N')
            low_conf_positions.append(i)
            continue
        # determine winner
        winner, winner_count = max(counts.items(), key=lambda x: x[1])
        freq = winner_count / total_non_gap
        if freq <= threshold:
            cons.append('N')
            low_conf_positions.append(i)
        else:
            cons.append(winner)
    # remove gaps from consensus
    consensus = ''.join(cons).replace('-', '')
    return consensus, low_conf_positions


def run():
    smk = globals().get("snakemake")
    if smk is None:
        raise SystemExit("This script is intended to be executed by Snakemake")

    aln = Path(smk.input["aln"]) if isinstance(smk.input, dict) else Path(smk.input[0])
    out_fa = Path(smk.output["fasta"]) if isinstance(smk.output, dict) else Path(smk.output[0])
    qc_path = Path(smk.output["qc"]) if isinstance(smk.output, dict) else None
    threshold = float(smk.params.get("threshold", 0.5))

    records = read_fasta(aln)
    cons_seq, low_conf = consensus_wta_with_conflict(records, threshold=threshold)

    out_fa.parent.mkdir(parents=True, exist_ok=True)
    sample = out_fa.stem.replace('.consensus', '')
    with open(out_fa, 'w') as out:
        out.write(f">{sample}\n")
        out.write(cons_seq + "\n")

    qc = {
        "input": str(aln),
        "output": str(out_fa),
        "threshold": threshold,
        "length": len(cons_seq),
        "n_aligned": len(records),
        "low_conf_positions": low_conf,
    }
    if qc_path:
        qc_path.parent.mkdir(parents=True, exist_ok=True)
        with open(qc_path, 'w') as f:
            json.dump(qc, f, indent=2)


if __name__ == "__main__":
    run()
