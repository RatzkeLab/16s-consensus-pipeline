# 16S Consensus Pipeline

A Snakemake pipeline for creating 16S consensus sequences from Oxford Nanopore (ONT) data, with support for multiple copy numbers per strain.

## Setup

### 1. Install Conda/Mamba
```bash
# Install Miniforge (includes mamba)
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```

### 2. Create Snakemake Environment
```bash
conda create -n snakemake_env -c conda-forge -c bioconda snakemake=9
conda activate snakemake_env
```

## Usage

### Basic Usage (Demo Data)
```bash
conda activate snakemake_env
snakemake -c 4 --use-conda
```

This runs the pipeline on demo data in `demo/merged/` with 4 cores.

### Custom Data
```bash
# 1. Create a config file (copy from demo/default_config.yaml)
cp demo/default_config.yaml my_config.yaml

# 2. Edit paths and parameters
nano my_config.yaml

# 3. Run pipeline
snakemake -c 8 --use-conda --configfile my_config.yaml
```

### Configuration Options

Key parameters in `demo/default_config.yaml`:

```yaml
input_dir: "demo/merged"              # Input FASTQ directory
output_dir: "demo/demo_output"        # Output directory
run_name: "demo"                      # Database filename prefix

min_reads_initial: 10                 # Min reads before filtering
min_reads_filtered: 5                 # Min reads after filtering

filter:
  min_avg_qscore: 10                  # NanoFilt quality threshold
  min_length: 1000                    # Minimum read length
  max_length: 2000                    # Maximum read length

subsample_n: 150                      # Reads per sample for alignment

consensus:
  min_consensus_proportion: 0.6       # 60% threshold for consensus
```

## Pipeline Steps

1. **Count Reads** - Count reads in each FASTQ file
2. **Initial QC Checkpoint** - Filter samples with ≥ min_reads_initial
3. **Filter Reads** - Apply NanoFilt quality/length filtering
4. **Filtered QC Checkpoint** - Filter samples with ≥ min_reads_filtered
5. **Subsample** - Randomly select up to subsample_n reads
6. **Align** - Create multiple sequence alignment with MAFFT
7. **Consensus** - Generate naive consensus (winner-takes-all, 60% threshold)
8. **Pool** - Concatenate all consensus sequences into single database FASTA
9. **Summary** - Generate pipeline report

## Outputs

```
output_dir/
  naive_db.fasta              # Naive consensus database (winner-takes-all)
  multi_db.fasta              # Multi-consensus database (with clustering)
  pipeline_summary.md         # QC report
  filtered/                   # Filtered FASTQs
  subsampled/                 # Subsampled FASTQs
  alignment/                  # Per-sample MSAs
  naive_consensus/            # Per-sample naive consensus FASTAs
  multi_consensus/            # Per-sample multi-consensus FASTAs (may have A, B, C variants)
  checks/                     # QC checkpoint files
  logs/                       # Per-sample logs
```

## Advanced Usage

**Visualize cluster separation (optional):**
```bash
conda activate snakemake_env
snakemake visualize -c 4 --use-conda
# Creates cluster_alignments/ with per-cluster alignment files
```

**Clean outputs and re-run:**
```bash
snakemake --delete-all-output
snakemake -c 4 --use-conda
```

**View summary report:**
```bash
cat output_dir/pipeline_summary.md
```

## Troubleshooting

**Dry-run to check workflow:**
```bash
snakemake -n
```

**Clean outputs and re-run:**
```bash
snakemake --delete-all-output
snakemake -c 4 --use-conda
```

**Check which samples failed QC:**
```bash
cat output_dir/pipeline_summary.md
```
