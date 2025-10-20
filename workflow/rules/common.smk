from pathlib import Path

# Base directories and config
OUT_DIR        = Path(config.get("outdir", "output"))
INPUT_DIR      = Path(config.get("fastq_input", "demo/merged"))

# Structured output directories
FILTER_DIR     = OUT_DIR / "filtered"
SAMPLE_DIR     = OUT_DIR / "sampled"
ALIGNMENT_DIR  = OUT_DIR / "alignment"
CONSENSUS_DIR  = OUT_DIR / "consensus"
QC_DIR         = OUT_DIR / "qc"
LOG_DIR        = OUT_DIR / "logs"

# Log subdirs
LOG_FILTER     = LOG_DIR / "filter"
LOG_ALIGNMENT  = LOG_DIR / "alignment"
LOG_SAMPLE     = LOG_DIR / "sample"
LOG_CONSENSUS  = LOG_DIR / "consensus"

# Final database file
CONSENSUS_DATABASE = OUT_DIR / config.get("database_filename", "db.fasta")


def load_bundle_names(input_dir: Path):
    """Return list of sample names from FASTQ files in input_dir.

    Supports both .fastq and .fastq.gz; takes unique stems.
    """
    names = []
    for pattern in ("*.fastq", "*.fastq.gz"):
        for p in sorted(input_dir.glob(pattern)):
            stem = p.name
            if stem.endswith(".fastq"):
                stem = p.stem
            elif stem.endswith(".fastq.gz"):
                stem = p.name[:-9]
            if stem:
                names.append(stem)
    names = sorted(set(names))
    if not names:
        raise ValueError(f"No *.fastq or *.fastq.gz files found in {input_dir}")
    return names

BUNDLE_NAMES = load_bundle_names(INPUT_DIR)


# Checkpoint helpers
def get_passing_samples(wildcards):
    """Load passing samples determined by the checkpoint.

    Reads names from qc/pass/passing_samples.txt.
    """
    pass_file = QC_DIR / "pass" / "passing_samples.txt"
    with open(pass_file) as fh:
        return [l.strip() for l in fh if l.strip()]


def passed_report_paths(wildcards):
    """Return per-sample report paths for samples that passed the checkpoint.

    This function is intended for use in the `input` section of rules (Snakemake 9 linting).
    """
    samples = get_passing_samples(wildcards)
    return [str(QC_DIR / "reports" / f"{s}.html") for s in samples]


def passed_alignment_qc(wildcards):
    samples = get_passing_samples(wildcards)
    return [str(QC_DIR / "alignment" / f"{s}.json") for s in samples]


def passed_consensus_qc(wildcards):
    samples = get_passing_samples(wildcards)
    return [str(QC_DIR / "consensus" / f"{s}.json") for s in samples]

