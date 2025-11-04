# Developer Guide

## Pipeline Architecture

### Overview
This is a Snakemake 9 pipeline that processes samples **independently** until the final aggregation step, maximizing parallelization. Each sample flows through checkpoints that filter out low-quality data before expensive operations.

### Directory Structure
workflow/
├── Snakefile # Main entry point, includes all rules
├── rules/
│ ├── common.smk # All paths, variables, helper functions
│ ├── 1_preprocessing.smk
│ ├── 2_align.smk
│ ├── 3_naive_consensus.smk
│ ├── 4_clustering.smk
│ ├── 5_cluster_consensus.smk
│ └── 6_summary.smk
├── scripts/
│ ├── checkpoint_helpers.py # Checkpoint aggregation functions
│ ├── common_helpers.py # Shared utility functions
│ ├── count_reads.py
│ ├── detect_clusters.py
│ ├── naive_consensus.py
│ ├── pairwise_distances.py
│ └── split_by_cluster.py
└── envs/ # Conda environment files

### Sample Flow
                             ┌→ naive_consensus ─┐
                             │                   ↓
                      align ─┤              pool_naive_db
                             │                   ↓
                             └→ detect_clusters  summary
                                     ↓
                                checkpoint
                                     ↓
                              split_by_cluster
                                     ↓
                               align_cluster
                                     ↓
                              cluster_consensus
                                     ↓
                               pool_multi_db


## Key Design Patterns

### 1. Checkpoints
Two critical checkpoints filter samples before expensive operations:
- **`check_min_reads_initial`** - before filtering
- **`check_min_reads_filtered`** - after filtering

These write `passing_samples.txt` files that downstream rules read via helper functions in `checkpoint_helpers.py`.

### 2. Parallel Consensus Paths
The pipeline generates two consensus databases simultaneously:
- **Naive path**: Simple plurality consensus (faster, always runs)
- **Multi path**: Cluster detection + per-cluster consensus (detects variants)

### 3. Single Source of Truth
All paths, variables, and parameters are declared in `common.smk` and `default_config.yaml`. Never hardcode paths in individual rules.

## Adding a New Rule

1. **Choose the right file** in `workflow/rules/`
2. **Add to common.smk** if you need new paths/variables
3. **Create conda env** in `workflow/envs/` if new dependencies needed
4. **Write minimal rule** - prefer external tools over custom scripts
5. **Update `rule all`** if it's a final output
6. **Test** with `snakemake -c 4 --use-conda -n` (dry-run)

Example:
````python
// filepath: workflow/rules/1_preprocessing.smk
rule example_rule:
    input:
        RESULTS_DIR / "{sample}" / "input.txt"
    output:
        RESULTS_DIR / "{sample}" / "output.txt"
    conda:
        "../envs/example.yaml"
    shell:
        "tool --input {input} --output {output}"
````

### Modifying Existing Rules
Adding Parameters
Add to demo/default_config.yaml
Access via config["your_param"]
Document in README.md
Changing Checkpoints
If you modify checkpoint logic, update the corresponding helper function in checkpoint_helpers.py to match.

### Testing
```bash
# Dry-run to check DAG
snakemake -c 4 --use-conda -n

# Run on demo data
conda activate snakemake_env
snakemake -c 4 --use-conda

# Clean and re-run
snakemake --delete-all-output
snakemake -c 4 --use-conda
```

### Common Patterns
Getting passing samples after a checkpoint
def get_samples_after_checkpoint(wildcards):
    checkpoint_output = checkpoints.checkpoint_name.get(**wildcards).output[0]
    return read_passing_samples(checkpoint_output.parent / "passing_samples.txt")

Aggregating per-sample outputs
input:
    expand(RESULTS_DIR / "{sample}" / "file.txt", sample=get_passing_samples)

Conditional execution
Use checkpoints to skip samples - don't use conditional logic in rules.

Snakemake 9 Specific
No ancient() - removed in v8+
Use Path objects from pathlib
Checkpoints must use checkpoints.name.get(**wildcards)
Avoid lambdas in inputs - use functions in common.smk
Quick Reference
Add software: Create/edit files in envs
Add config option: Edit default_config.yaml
Add output path: Declare in common.smk
Skip samples: Use checkpoints, not conditionals
Parallelize: Keep per-sample until aggregation


This gives developers the mental model and patterns they need without duplicating user-facing info. Let me know if you want any sections expanded or condensed!This gives developers the mental model and patterns they need without duplicating user-facing info. Let me know if you want any sections expanded or condensed!