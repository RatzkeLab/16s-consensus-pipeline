"""
Alignment, consensus, and consolidation rules.
"""

rule align_sample:
	input:
		fasta = str(SAMPLE_DIR / "{sample}.sampled.fasta")
	output:
		aln = temp(str(ALIGNMENT_DIR / "{sample}.aln.fasta")),
		qc  = str(QC_DIR / "alignment" / "{sample}.json")
	log:
		str(LOG_ALIGNMENT / "{sample}.log")
	conda: "../envs/align.yaml"
	shell:
		r"""
		set -euo pipefail
		mkdir -p "$(dirname {output.aln})" "$(dirname {output.qc})" "$(dirname {log})"
		if command -v mafft >/dev/null 2>&1; then
			mafft --auto {input.fasta} > {output.aln} 2> {log}
			echo '{{"tool":"mafft"}}' > {output.qc}
		else
			# fallback: copy as single-sequence alignment
			cp {input.fasta} {output.aln}
			echo '{{"tool":"none"}}' > {output.qc}
		fi
		"""


rule consensus_from_alignment:
	input:
		aln = str(ALIGNMENT_DIR / "{sample}.aln.fasta")
	output:
		fasta = str(CONSENSUS_DIR / "{sample}.consensus.fasta"),
		qc    = str(QC_DIR / "consensus" / "{sample}.json")
	params:
		threshold = 0.5
	conda: "../envs/qc.yaml"
	script:
		"../scripts/consensus.py"


rule build_database:
	input:
		expand(str(CONSENSUS_DIR / "{sample}.consensus.fasta"), sample=BUNDLE_NAMES)
	output:
		db = str(CONSENSUS_DATABASE)
	conda: "../envs/qc.yaml"
	shell:
		r"""
		mkdir -p "$(dirname {output.db})"
		cat {input} > {output.db}
		"""


rule report_alignment_summary:
	input:
		expand(str(QC_DIR / "alignment" / "{sample}.json"), sample=BUNDLE_NAMES)
	output:
		tsv = str(QC_DIR / "reports" / "alignment_summary.tsv")
	conda: "../envs/qc.yaml"
	script:
		"../scripts/report_alignment_summary.py"


rule report_consensus_summary:
	input:
		expand(str(QC_DIR / "consensus" / "{sample}.json"), sample=BUNDLE_NAMES)
	output:
		tsv = str(QC_DIR / "reports" / "consensus_summary.tsv")
	conda: "../envs/qc.yaml"
	script:
		"../scripts/report_consensus_summary.py"


rule per_sample_report:
	input:
		filter_qc   = str(QC_DIR / "filter" / "{sample}.json"),
		sample_qc   = str(QC_DIR / "sample" / "{sample}.json"),
		align_qc    = str(QC_DIR / "alignment" / "{sample}.json"),
		consensus_qc= str(QC_DIR / "consensus" / "{sample}.json")
	output:
		html = str(QC_DIR / "reports" / "{sample}.html")
	conda: "../envs/qc.yaml"
	script:
		"../scripts/per_sample_report.py"

