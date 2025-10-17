
from pathlib import Path

INPUT_DIR       = Path(config.get("fastq_input", "demo/fastq_pass_partial"))
FILTER_DIR      = OUT_DIR / "filtered"
LOG_FILTER      = OUT_DIR / "logs" / "filter"
ALIGNMENT_DIR   = OUT_DIR / "alignment"
LOG_ALIGNMENT   = OUT_DIR / "logs" / "alignment"
SPLIT_DIR      = OUT_DIR / "split"
LOG_SPLIT       = OUT_DIR / "logs" / "split"


CONSENSUS_DATABASE = OUT_DIR / config.get("database_filename", "db.fasta")