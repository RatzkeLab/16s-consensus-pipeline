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
conda create -n snakemake_env -c conda-forge -c bioconda snakemake
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
snakemake -c 8 --use-conda --configfile path/to/your/config/my_config.yaml
```

### Configuration Options

Key parameters that should be modified in your config:

```yaml
input_dir: "path/to/your/fastq/folder"              # Input FASTQ directory
output_dir: "path/to/wherever/you/want/your/output/to/be/saved"        # Output directory
run_name: "demo"                      # Database filename prefix

min_reads_initial: 10                 # Min reads before filtering
min_reads_filtered: 5                 # Min reads after filtering

filter:
  min_avg_qscore: 10                  # NanoFilt quality threshold
  min_length: 1000                    # Minimum read length
  max_length: 2000                    # Maximum read length

subsample_n: 150                      # Reads per sample for alignment

consensus:
  min_consensus_proportion: 0.6       # 60% threshold for adding a locus to the error profile
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

Outputs are generated along the way at each step, and can be used to sanity check the performance of the pipeline on your data.

The key final outputs are 
 - naive_db.fasta, which contains 1 consensus sequence per sample
 - multi_db.fasta, which contains multiple consensus sequences per sample, if multiple consensuses are detected.

## Advanced Usage

**Clean outputs and re-run:**
```bash
snakemake --delete-all-output
snakemake -c 4 --use-conda
```