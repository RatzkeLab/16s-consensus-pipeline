"""
Sampling rules and pass-samples checkpoint.
"""

checkpoint pass_samples:
	input:
		expand(str(QC_DIR / "filter" / "{sample}.json"), sample=BUNDLE_NAMES)
	output:
		pass_list = str(QC_DIR / "pass" / "passing_samples.txt")
	params:
		min_reads = 10
	conda: "../envs/qc.yaml"
	script:
		"../scripts/determine_passing_samples.py"


rule sample_reads:
	input:
		fastq = str(FILTER_DIR / "{sample}.filtered.fastq.gz")
	output:
		fasta = temp(str(SAMPLE_DIR / "{sample}.sampled.fasta")),
		qc    = str(QC_DIR / "sample" / "{sample}.json")
	params:
		n = int(config["alignment"].get("n_sequences", 150)),
		method = config["alignment"].get("sample_method", "random")
	log:
		str(LOG_SAMPLE / "{sample}.log")
	conda: "../envs/qc.yaml"
	script:
		"../scripts/sample_reads.py"


def samples_after_checkpoint(wildcards):
	return get_passing_samples(wildcards)


