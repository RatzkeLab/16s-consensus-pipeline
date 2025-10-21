import sys, collections, csv
from pathlib import Path

aln_path = "{input.aln}"
fa_out   = "{output.fa}"
csv_out  = "{output.csv}"
sample   = "{wildcards.sample}"

# read FASTA
seqs = []
with open(aln_path) as f:
    name=None; buf=[]
    for line in f:
        line=line.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None:
                seqs.append(("".join(buf)).upper())
            name = line[1:]; buf=[]
        else:
            buf.append(line)
    if name is not None:
        seqs.append(("".join(buf)).upper())

if not seqs:
    print(f"ERROR: empty alignment {aln_path}", file=sys.stderr)
    sys.exit(1)

L = len(seqs[0])
if any(len(s)!=L for s in seqs):
    print("ERROR: alignment not equal-length", file=sys.stderr); sys.exit(1)

valid = set("ACGT")
cons = []
low_rows = []  # columns where top < 0.5

for i in range(L):
    col = [s[i] for s in seqs if s[i] in valid]
    n = len(col)
    if n == 0:
        cons.append("N"); continue
    cnt = collections.Counter(col)
    top_base, top_count = cnt.most_common(1)[0]
    prop = top_count / n
    # plurality by definition (no >=0.5 requirement); ties handled below
    # if tie for top, emit N
    # detect tie
    mult = cnt.most_common()
    if len(mult) > 1 and mult[0][1] == mult[1][1]:
        cons.append("N")
        low_rows.append((i+1, "N", mult[0][1], n, mult[0][1]/n, "tie"))
        continue
    cons.append(top_base)
    if prop < 0.5:
        low_rows.append((i+1, top_base, top_count, n, prop, ""))

# write FASTA (strip gaps)
seq = "".join(cons).replace("-", "")
with open(fa_out, "w") as g:
    g.write(f">{sample}\n")
    for j in range(0, len(seq), 80):
        g.write(seq[j:j+80] + "\n")

# write CSV of <50% consensus sites (and ties as flagged)
with open(csv_out, "w", newline="") as g:
    w = csv.writer(g)
    w.writerow(["position","consensus","top_count","column_non_gap","prop_top","note"])
    for r in low_rows:
        w.writerow(r)