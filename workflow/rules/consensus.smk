"""
Naive consensus sequence generation from alignments.
"""

rule naive_consensus:
    """
    Generate naive consensus sequence using winner-takes-all approach.
    Positions below threshold are marked with winner base and variants recorded.
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        fasta=NAIVE_CONSENSUS_DIR / "{sample}.fasta",
        variants=NAIVE_CONSENSUS_DIR / "{sample}_variants.tsv"
    params:
        sample="{sample}",
        min_prop=NAIVE_CONSENSUS_MIN_PROP
    log:
        LOG_DIR / "naive_consensus" / "{sample}.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p "$(dirname {output.fasta})"
        
        python workflow/scripts/naive_consensus.py \
          {input.alignment} \
          {output.fasta} \
          {output.variants} \
          {params.sample} \
          {params.min_prop} \
          2> {log}
        """