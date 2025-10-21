"""
Consensus sequence generation from alignments.
"""

rule consensus:
    """
    Generate naive consensus sequence using winner-takes-all approach.
    Positions below threshold are marked as N and variants recorded.
    """
    input:
        alignment=ALIGNMENT_DIR / "{sample}.fasta"
    output:
        fasta=CONSENSUS_DIR / "{sample}.fasta",
        variants=CONSENSUS_DIR / "{sample}_variants.tsv"
    params:
        sample="{sample}",
        min_prop=CONSENSUS_MIN_PROP
    log:
        LOG_DIR / "consensus" / "{sample}.log"
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