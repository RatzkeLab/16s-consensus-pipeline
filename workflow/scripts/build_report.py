import os, glob, html, pathlib, csv

title = "{TITLE}"
outdir = "{OUT}"

# collect inputs by suffix
cons = sorted([p for p in {input} if p.endswith(".fasta") and "/consensus/" in p])
stats= sorted([p for p in {input} if p.endswith(".csv") and "/consensus_stats/" in p])
pre  = sorted([p for p in {input} if p.endswith("_pre.tsv")])
post = sorted([p for p in {input} if p.endswith("_post.tsv")])

def read_tsv_one(p):
    row = list(csv.DictReader(open(p, newline=''), delimiter="\t"))[0]
    return row

def fasta_len(p):
    L=0
    with open(p) as f:
        for line in f:
            if not line.startswith(">"):
                L += len(line.strip())
    return L

# maps for counts
pre_map = {pathlib.Path(p).stem.replace("_pre",""): read_tsv_one(p) for p in pre}
post_map= {pathlib.Path(p).stem.replace("_post",""): read_tsv_one(p) for p in post}
low_map = {}
for p in stats:
    s = pathlib.Path(p).stem
    # count rows
    n = sum(1 for _ in open(p)) - 1 if os.path.getsize(p) > 0 else 0
    low_map[s] = max(n, 0)

rows=[]
for p in cons:
    s = pathlib.Path(p).stem
    rows.append((
        s,
        pre_map.get(s,{}).get("reads","NA"),
        pre_map.get(s,{}).get("status","NA"),
        post_map.get(s,{}).get("reads","NA"),
        post_map.get(s,{}).get("status","NA"),
        fasta_len(p),
        low_map.get(s, 0)
    ))
rows.sort()

print("<!doctype html><meta charset='utf-8'>")
print(f"<title>{html.escape(title)}</title>")
print(f"<h1>{html.escape(title)}</h1>")
print("<p>Consensus generated only for samples that passed both read thresholds.</p>")
print("<table border='1' cellspacing='0' cellpadding='4'>")
print("<tr><th>Sample</th><th>Reads pre</th><th>Status</th><th>Reads post</th><th>Status</th><th>Consensus length</th><th>#Sites &lt;50%</th></tr>")
for r in rows:
    s, pre_n, pre_st, post_n, post_st, clen, lowc = r
    print(f"<tr><td>{html.escape(s)}</td><td>{pre_n}</td><td>{pre_st}</td><td>{post_n}</td><td>{post_st}</td><td>{clen}</td><td>{lowc}</td></tr>")
print("</table>")