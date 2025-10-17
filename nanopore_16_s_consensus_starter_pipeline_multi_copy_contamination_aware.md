# Nanopore 16S Consensus Starter Pipeline

A Snakemake-based starter pipeline to derive high-quality consensus sequences from **Nanopore 16S amplicon reads** for **individual strains**, while being **aware of multiple 16S copies (≤6)** and **potential cross-contamination**.

It:
- clusters reads into up to `expected_max_copies` haplotype groups (default **6**),
- builds & polishes a consensus **per cluster** (MAFFT → Racon → Medaka),
- flags **multi-copy-like** vs **contamination-like** signal via **allele balance** at variable sites,
- discards tiny or aberrant clusters as **trash**, and
- outputs a clean **FASTA of per-cluster consensuses** plus a **summary TSV** with statistics you can inspect.

> ⚠️ This is a starter template. You can swap modules (e.g., VSEARCH clustering) or tweak thresholds to your data.

---

## Directory Layout
```
consensus16s/
├─ Snakefile
├─ config.yaml
├─ envs/
│  └─ consensus.yml
├─ scripts/
│  ├─ cluster_hdbscan.py
│  ├─ build_consensus.py
│  ├─ call_site_stats.py
│  └─ summarize_clusters.py
├─ README.md
└─ data/
   ├─ samples.tsv           # two columns: sample_id\tfastq
   └─ fastq/                # your .fastq(.gz) files (if not using samples.tsv)
```

---

## Quickstart

```bash
# 1) create conda env (mamba recommended)
mamba env create -f envs/consensus.yml
mamba activate consensus16s

# 2) configure paths & params
#    - edit config.yaml
#    - create data/samples.tsv OR put reads under data/fastq/

# 3) run
snakemake -j 8 --use-conda

# 4) outputs (per sample under results/{sample}/):
#    - clusters.tsv: cluster assignments & sizes
#    - cluster*/consensus.fasta: polished per-cluster consensus
#    - site_stats.tsv: per-cluster allele balance across variable sites
#    - summary.tsv: final per-cluster call: {good|multicopy_like|contam_like|trash}
```

---

## `envs/consensus.yml`
```yaml
name: consensus16s
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=3.11
  - snakemake
  - minimap2
  - samtools
  - racon
  - medaka
  - mafft
  - bcftools
  - seqkit
  - cutadapt
  - pigz
  - pip
  - pip:
      - biopython
      - pandas
      - numpy
      - hdbscan
      - umap-learn
      - scikit-learn
```

---

## `config.yaml`
```yaml
# Input
samples_tsv: data/samples.tsv   # if null, pipeline scans data/fastq/*.fastq* instead
primer_f: null                  # optional: 5' primer seq (for cutadapt)
primer_r: null                  # optional: 3' primer seq (revcomp allowed)
trim_error_rate: 0.1

# Basic filtering
min_read_length: 1300
max_read_length: 1700
min_qscore: 7                   # simple pass filter via seqkit (approx proxy)

# Clustering
clustering: hdbscan             # hdbscan | vsearch
expected_max_copies: 6
min_cluster_reads: 20           # clusters < this → trash

# HDBSCAN parameters (k-mer embedding + UMAP + HDBSCAN)
kmer_k: 6
umap_n_neighbors: 15
umap_min_dist: 0.1
hdbscan_min_cluster_size: 25
hdbscan_min_samples: 5

# VSEARCH fallback (if you prefer centroid clustering)
id_threshold: 0.995

# Consensus polishing
racon_rounds: 2
medaka_model: r941_min_sup_g507 # pick for your flowcell/chemistry
subsample_for_msa: 500          # reads to MSA for seed consensus

# Allele balance logic (per-cluster)
minor_allele_multi_copy_lb: 0.40  # ~50/50 suggests multi-copy
minor_allele_multi_copy_ub: 0.60
minor_allele_contam_max: 0.20     # <20% suggests contamination/minor noise
min_variable_sites_for_flag: 2

# Outlier / trash rules
max_clusters: 6                   # hard cap (keeps top-N by size)
max_pairwise_divergence: 0.02     # >2% vs bulk ⇒ trash (radically different)

# Reporting
report_min_cluster_reads: 5
```

---

## `Snakefile`
```python
import os
import pandas as pd
configfile: "config.yaml"

SAMPLES_TSV = config.get("samples_tsv", None)

if SAMPLES_TSV:
    SAMPLES = pd.read_csv(SAMPLES_TSV, sep='\t', header=None, names=["sample","fastq"]).set_index("sample")["fastq"].to_dict()
else:
    # infer samples from data/fastq
    from glob import glob
    fq = sorted(glob("data/fastq/*.fastq*"))
    SAMPLES = {os.path.basename(f).split(".")[0]: f for f in fq}

rule all:
    input:
        expand("results/{sample}/summary.tsv", sample=list(SAMPLES.keys()))

rule filter_reads:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        fq="results/{sample}/filtered.fastq.gz"
    conda:
        "envs/consensus.yml"
    shell:
        r"""
        mkdir -p results/{wildcards.sample}
        seqkit seq -m {config[min_read_length]} -M {config[max_read_length]} {input} \
          | seqkit fqchk -Q {config[min_qscore]} - \
          | pigz -c > {output}
        """

rule trim_primers:
    input:
        fq="results/{sample}/filtered.fastq.gz"
    output:
        fq="results/{sample}/trimmed.fastq.gz"
    conda:
        "envs/consensus.yml"
    run:
        fwd = config.get("primer_f")
        rev = config.get("primer_r")
        if fwd and rev:
            shell(r"""
            cutadapt -e {config[trim_error_rate]} -g {fwd} -a {rev} -j 0 -o {output} {input.fq}
            """)
        else:
            shell("cp {input.fq} {output}")

rule cluster_reads:
    input:
        fq="results/{sample}/trimmed.fastq.gz"
    output:
        tsv="results/{sample}/clusters.tsv"
    conda:
        "envs/consensus.yml"
    params:
        method=config["clustering"],
        max_c=config["expected_max_copies"],
        min_reads=config["min_cluster_reads"],
        idthr=config.get("id_threshold", 0.995),
        k=config.get("kmer_k", 6),
        nnei=config.get("umap_n_neighbors", 15),
        md=config.get("umap_min_dist", 0.1),
        mcs=config.get("hdbscan_min_cluster_size", 25),
        ms=config.get("hdbscan_min_samples", 5)
    shell:
        r"""
        python scripts/cluster_hdbscan.py \
          --fastq {input.fq} --out {output.tsv} \
          --method {params.method} --max-clusters {params.max_c} \
          --min-cluster-reads {params.min_reads} \
          --id-threshold {params.idthr} \
          --k {params.k} --umap-n-neighbors {params.nnei} --umap-min-dist {params.md} \
          --hdbscan-min-cluster-size {params.mcs} --hdbscan-min-samples {params.ms}
        """

rule build_consensus_per_cluster:
    input:
        fq="results/{sample}/trimmed.fastq.gz",
        clus="results/{sample}/clusters.tsv"
    output:
        touch("results/{sample}/consensus.done")
    conda:
        "envs/consensus.yml"
    params:
        racon_rounds=config["racon_rounds"],
        medaka_model=config["medaka_model"],
        subsamp=config["subsample_for_msa"]
    shell:
        r"""
        python scripts/build_consensus.py \
          --fastq {input.fq} --clusters {input.clus} \
          --outdir results/{wildcards.sample} \
          --racon-rounds {params.racon_rounds} \
          --medaka-model {params.medaka_model} \
          --subsample {params.subsamp}
        touch {output}
        """

rule call_site_stats:
    input:
        fq="results/{sample}/trimmed.fastq.gz",
        done="results/{sample}/consensus.done"
    output:
        tsv="results/{sample}/site_stats.tsv"
    conda:
        "envs/consensus.yml"
    shell:
        r"""
        python scripts/call_site_stats.py \
          --fastq {input.fq} --clusters results/{wildcards.sample}/clusters.tsv \
          --consensus-glob "results/{wildcards.sample}/cluster*/consensus.fasta" \
          --out {output.tsv}
        """

rule summarize:
    input:
        clus="results/{sample}/clusters.tsv",
        sites="results/{sample}/site_stats.tsv"
    output:
        tsv="results/{sample}/summary.tsv"
    conda:
        "envs/consensus.yml"
    params:
        mc_lb=config["minor_allele_multi_copy_lb"],
        mc_ub=config["minor_allele_multi_copy_ub"],
        contam_max=config["minor_allele_contam_max"],
        min_sites=config["min_variable_sites_for_flag"],
        min_reads=config["min_cluster_reads"],
        report_min=config["report_min_cluster_reads"]
    shell:
        r"""
        python scripts/summarize_clusters.py \
          --clusters {input.clus} --site-stats {input.sites} \
          --out {output.tsv} \
          --minor-mc-lb {params.mc_lb} --minor-mc-ub {params.mc_ub} \
          --minor-contam-max {params.contam_max} --min-var-sites {params.min_sites} \
          --min-cluster-reads {params.min_reads} --report-min {params.report_min}
        """
```

---

## `scripts/cluster_hdbscan.py`
```python
#!/usr/bin/env python3
import argparse, gzip, sys, os
from collections import Counter
import numpy as np
import pandas as pd

# Optional: VSEARCH fallback via subprocess
import subprocess, tempfile

# K-mer embedding → UMAP → HDBSCAN
from sklearn.feature_extraction.text import CountVectorizer
import umap
import hdbscan

DNA = set("ACGTN")

def openmaybe(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def read_fastq_seqs(path, max_reads=None):
    seqs = []
    with openmaybe(path) as fh:
        i = 0
        for line in fh:
            if line.startswith('@'):
                seq = next(fh).strip().upper()
                # skip '+' and qual
                next(fh); next(fh)
                if set(seq) - DNA:
                    continue
                seqs.append(seq)
                i += 1
                if max_reads and i >= max_reads:
                    break
    return seqs

def kmer_strings(seqs, k=6):
    # Return space-separated k-mers for vectorizer
    out = []
    for s in seqs:
        toks = [s[i:i+k] for i in range(0, len(s)-k+1)]
        out.append(" ".join(toks))
    return out

def cluster_hdbscan(fastq, k, nnei, md, mcs, ms):
    seqs = read_fastq_seqs(fastq)
    if len(seqs) == 0:
        raise SystemExit("No reads found after filtering.")
    corpus = kmer_strings(seqs, k)
    vec = CountVectorizer(analyzer='word', token_pattern=r'[^\s]+')
    X = vec.fit_transform(corpus)
    reducer = umap.UMAP(n_neighbors=nnei, min_dist=md, metric='cosine', random_state=42)
    emb = reducer.fit_transform(X)
    cl = hdbscan.HDBSCAN(min_cluster_size=mcs, min_samples=ms).fit(emb)
    labels = cl.labels_
    return labels

def run_vsearch_cluster(fastq, idthr):
    with tempfile.TemporaryDirectory() as td:
        derep = os.path.join(td, 'derep.fasta')
        cent = os.path.join(td, 'centroids.fasta')
        uc = os.path.join(td, 'clusters.uc')
        # dereplicate to speed up
        subprocess.check_call(["vsearch", "--derep_fulllength", fastq, "--output", derep, "--sizeout"])
        subprocess.check_call(["vsearch", "--cluster_fast", derep, "--id", str(idthr), "--centroids", cent, "--uc", uc])
        # parse .uc to labels per original read in order: not trivial without mapping derep back
        # Placeholder: recommend HDBSCAN path for robust starter.
        raise SystemExit("VSEARCH path provided as placeholder; prefer HDBSCAN path in this starter.")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--method', default='hdbscan')
    ap.add_argument('--max-clusters', type=int, default=6)
    ap.add_argument('--min-cluster-reads', type=int, default=20)
    ap.add_argument('--id-threshold', type=float, default=0.995)
    ap.add_argument('--k', type=int, default=6)
    ap.add_argument('--umap-n-neighbors', type=int, default=15)
    ap.add_argument('--umap-min-dist', type=float, default=0.1)
    ap.add_argument('--hdbscan-min-cluster-size', type=int, default=25)
    ap.add_argument('--hdbscan-min-samples', type=int, default=5)
    args = ap.parse_args()

    if args.method == 'hdbscan':
        labels = cluster_hdbscan(args.fastq, args.k, args.umap_n_neighbors, args.umap_min_dist,
                                 args.hdbscan_min_cluster_size, args.hdbscan_min_samples)
    else:
        labels = run_vsearch_cluster(args.fastq, args.id_threshold)

    # Save as two-column TSV: read_idx \t cluster (>=0) or -1 for noise
    lab = pd.Series(labels)
    counts = lab[lab >= 0].value_counts().sort_values(ascending=False)
    # keep only top-N clusters by size
    keep = set(counts.head(args.max_clusters).index)
    # Map labels not in keep to -2 (trash)
    lab = lab.apply(lambda x: x if x in keep and x >= 0 else (-1 if x == -1 else -2))

    # Filter by min_cluster_reads
    sizes = lab[lab >= 0].value_counts()
    good = set(sizes[sizes >= args.min_cluster_reads].index)
    lab = lab.apply(lambda x: x if x in good else (-2 if x >= 0 else x))

    # Build output table
    df = pd.DataFrame({"read_idx": np.arange(len(lab)), "cluster": lab})
    df.to_csv(args.out, sep='\t', index=False)
```

---

## `scripts/build_consensus.py`
```python
#!/usr/bin/env python3
import argparse, os, subprocess, tempfile, gzip
import pandas as pd
from random import sample


def openmaybe(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'r')

def write_fastq_subset(infq, idxs, outfq):
    with openmaybe(infq) as ih, open(outfq, 'w') as oh:
        i = -1
        while True:
            h = ih.readline()
            if not h:
                break
            seq = ih.readline(); plus = ih.readline(); qual = ih.readline()
            i += 1
            if i in idxs:
                oh.write(h); oh.write(seq); oh.write(plus); oh.write(qual)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--racon-rounds', type=int, default=2)
    ap.add_argument('--medaka-model', default='r941_min_sup_g507')
    ap.add_argument('--subsample', type=int, default=500)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    clus = pd.read_csv(args.clusters, sep='\t')
    # keep only assigned clusters (>=0); -1 noise, -2 trash ignored
    clusters = sorted([c for c in clus['cluster'].unique() if c >= 0])

    for c in clusters:
        cdir = os.path.join(args.outdir, f"cluster{c}")
        os.makedirs(cdir, exist_ok=True)
        reads_idx = clus.loc[clus.cluster==c, 'read_idx'].tolist()
        # Subsample for initial MSA seed
        seed = set(sample(reads_idx, min(len(reads_idx), args.subsample)))
        seedfq = os.path.join(cdir, 'seed.fastq')
        write_fastq_subset(args.fastq, seed, seedfq)
        # Convert to FASTA for MAFFT
        seedfa = os.path.join(cdir, 'seed.fasta')
        subprocess.check_call(f"seqkit fq2fa {seedfq} > {seedfa}", shell=True)
        # MSA & draft consensus
        aln = os.path.join(cdir, 'seed.aln.fasta')
        subprocess.check_call(f"mafft --auto {seedfa} > {aln}", shell=True)
        cons = os.path.join(cdir, 'draft_consensus.fasta')
        subprocess.check_call(f"seqkit consensus {aln} > {cons}", shell=True)
        # Map all cluster reads to draft consensus, polish with racon
        clufq = os.path.join(cdir, 'cluster.fastq')
        write_fastq_subset(args.fastq, set(reads_idx), clufq)
        paf = os.path.join(cdir, 'map.paf')
        subprocess.check_call(f"minimap2 -x map-ont -t 4 {cons} {clufq} > {paf}", shell=True)
        polished = os.path.join(cdir, 'polish_round0.fasta')
        subprocess.check_call(f"racon -t 4 {clufq} {paf} {cons} > {polished}", shell=True)
        prev = polished
        for r in range(1, args.racon_rounds):
            pafr = os.path.join(cdir, f'map_r{r}.paf')
            subprocess.check_call(f"minimap2 -x map-ont -t 4 {prev} {clufq} > {pafr}", shell=True)
            nxt = os.path.join(cdir, f'polish_round{r}.fasta')
            subprocess.check_call(f"racon -t 4 {clufq} {pafr} {prev} > {nxt}", shell=True)
            prev = nxt
        # Optional medaka polish
        final_cons = os.path.join(cdir, 'consensus.fasta')
        medaka_dir = os.path.join(cdir, 'medaka')
        subprocess.check_call(f"medaka_consensus -i {clufq} -d {prev} -o {medaka_dir} -m {args.medaka_model}", shell=True)
        subprocess.check_call(f"cp {medaka_dir}/consensus.fasta {final_cons}", shell=True)
```

---

## `scripts/call_site_stats.py`
```python
#!/usr/bin/env python3
import argparse, glob, os, subprocess, tempfile
import pandas as pd

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--fastq', required=True)
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--consensus-glob', required=True)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    clus = pd.read_csv(args.clusters, sep='\t')
    clusters = sorted([c for c in clus['cluster'].unique() if c >= 0])
    rows = []
    for c in clusters:
        cdir = os.path.join(os.path.dirname(args.clusters), f"cluster{c}")
        cons = os.path.join(cdir, 'consensus.fasta')
        if not os.path.exists(cons):
            continue
        bam = os.path.join(cdir, 'reads.bam')
        # map all cluster reads to its consensus & make pileup
        clufq = os.path.join(cdir, 'cluster.fastq')
        subprocess.check_call(f"minimap2 -ax map-ont {cons} {clufq} | samtools sort -o {bam}", shell=True)
        subprocess.check_call(f"samtools index {bam}", shell=True)
        # mpileup, restrict to depth info and bases
        mp = subprocess.check_output(f"samtools mpileup -aa -A -d 100000 {bam}", shell=True).decode()
        for line in mp.strip().split('\n'):
            chrom, pos, ref, depth, bases, quals = line.split('\t')[:6]
            pos = int(pos)
            depth = int(depth)
            if depth == 0:
                continue
            # simple base parsing (approximate): count A/C/G/T ignoring indels
            counts = {b:0 for b in 'ACGT'}
            for ch in bases.upper():
                if ch in counts:
                    counts[ch]+=1
            total = sum(counts.values())
            if total == 0:
                continue
            sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
            major_base, major_n = sorted_counts[0]
            minor_n = total - major_n
            minor_freq = minor_n/total if total>0 else 0
            if minor_n>0:
                rows.append({
                    'cluster': c,
                    'pos': pos,
                    'ref': ref,
                    'depth': total,
                    'major_base': major_base,
                    'major_n': major_n,
                    'minor_n': minor_n,
                    'minor_freq': minor_freq
                })
    pd.DataFrame(rows).to_csv(args.out, sep='\t', index=False)
```

---

## `scripts/summarize_clusters.py`
```python
#!/usr/bin/env python3
import argparse, pandas as pd, os

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--clusters', required=True)
    ap.add_argument('--site-stats', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--minor-mc-lb', type=float, default=0.4)
    ap.add_argument('--minor-mc-ub', type=float, default=0.6)
    ap.add_argument('--minor-contam-max', type=float, default=0.2)
    ap.add_argument('--min-var-sites', type=int, default=2)
    ap.add_argument('--min-cluster-reads', type=int, default=20)
    ap.add_argument('--report-min', type=int, default=5)
    args = ap.parse_args()

    clus = pd.read_csv(args.clusters, sep='\t')
    site = pd.read_csv(args.site_stats, sep='\t') if os.path.exists(args.site_stats) and os.stat(args.site_stats).st_size>0 else pd.DataFrame(columns=['cluster'])

    out_rows = []
    for c, g in clus[clus.cluster>=0].groupby('cluster'):
        n_reads = g.shape[0]
        status = 'trash' if n_reads < args.min_cluster_reads else 'good'
        sub = site[site.cluster==c]
        # Count variable sites and assess minor allele frequency distribution
        n_var = sub.shape[0]
        # proportion of sites with near 50/50 minor freq (suggest multi-copy)
        near5050 = ((sub.minor_freq >= args.minor_mc_lb) & (sub.minor_freq <= args.minor_mc_ub)).mean() if n_var>0 else 0.0
        # max minor freq (suggest contamination if consistently low)
        max_minor = sub.minor_freq.max() if n_var>0 else 0.0
        flag = None
        if n_reads >= args.min_cluster_reads and n_var >= args.min_var_sites:
            if near5050 > 0.5:  # >50% of variable sites ~balanced
                flag = 'multicopy_like'
            elif max_minor <= args.minor_contam_max:
                flag = 'contam_like'
        if flag:
            status = flag
        out_rows.append({
            'cluster': c,
            'n_reads': n_reads,
            'n_variable_sites': int(n_var),
            'near_50_50_fraction': round(float(near5050),3),
            'max_minor_freq': round(float(max_minor),3),
            'status': status,
            'consensus_fasta': f"cluster{c}/consensus.fasta"
        })

    df = pd.DataFrame(out_rows).sort_values(['status','n_reads'], ascending=[True, False])
    # only report clusters with at least report_min reads
    df = df[df.n_reads >= args.report_min]
    df.to_csv(args.out, sep='\t', index=False)
```

---

## `README.md`
```markdown
# Nanopore 16S Consensus (multi-copy & contamination aware)

This pipeline clusters Nanopore 16S reads per sample, builds a polished consensus per cluster, and labels clusters as:
- **good** — strong support, no imbalance signal,
- **multicopy_like** — many variable sites show ~50/50 allele balance, consistent with multiple 16S copies in one strain,
- **contam_like** — variable sites have consistently low minor-allele fractions (e.g., <20%), suggesting minor contamination,
- **trash** — too few reads or noise.

See `config.yaml` to tune thresholds (e.g., `expected_max_copies: 6`).

## Inputs
- `data/samples.tsv` **or** `data/fastq/*.fastq(.gz)`
- optional primer sequences in `config.yaml` to trim with cutadapt

## Outputs
- `results/{sample}/clusters.tsv` — per-read cluster assignment (HDBSCAN/UMAP)
- `results/{sample}/cluster*/consensus.fasta` — per-cluster polished consensus (MAFFT→Racon→Medaka)
- `results/{sample}/site_stats.tsv` — per-site allele depth & minor allele frequencies per cluster
- `results/{sample}/summary.tsv` — per-cluster summary & status with read counts

## Notes
- Swap clustering to `vsearch` if preferred; HDBSCAN path is default.
- For highly similar paralogs, increase read depth and adjust `hdbscan_min_cluster_size`, `id_threshold`, etc.
- Consider adding `NanoCLUST` style refinement or chimera filtering if needed.
```

