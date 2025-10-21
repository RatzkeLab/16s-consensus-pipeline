"""
Multiple sequence alignment using MAFFT.
"""

rule align:
    """
    Create multiple sequence alignment per sample using MAFFT.
    """
    input:
        fastq=SUBSAMPLE_DIR / "{sample}.fastq"
    output:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    log:
        LOG_DIR / "align" / "{sample}.log"
    conda:
        "../envs/align.yaml"
    threads: 4
    shell:
        """
        mkdir -p "$(dirname {output.alignment})"
        
        # Convert FASTQ to FASTA and align with MAFFT
        awk 'NR%4==1{{print ">"substr($0,2)}} NR%4==2{{print}}' {input.fastq} \
          | mafft --thread {threads} --auto - \
          > {output.alignment} 2> {log}
        """
