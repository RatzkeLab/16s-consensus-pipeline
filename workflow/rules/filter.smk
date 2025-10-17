#implicit import of common.smk for helpers and shared globals

rule filter:
    input:
        str(INPUT_DIR / "{read_bundle}.fastq.gz")
    output:
        temp(str(FILTER_DIR / "{read_bundle}.filtered.fastq.gz"))
    params:
        min_len = config["nanofilt"]["min_len"],
        min_q   = config["nanofilt"]["min_q"],
        head    = config["nanofilt"].get("headcrop", 0),
        tail    = config["nanofilt"].get("tailcrop", 0),
        max_len_str = (f'--maxlength {config["nanofilt"]["max_len"]}') if config["nanofilt"].get("max_len", -1) > 0 else '',
    log:
        str(LOG_FILTER / "{read_bundle}.log")
    conda: "../envs/nanofilt.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname "{output}")" "$(dirname "{log}")"
        gzip -cd {input} \
        | NanoFilt --length {params.min_len} --quality {params.min_q} \
                   --headcrop {params.head} --tailcrop {params.tail} \
                   {params.max_len_str} \
                   2> "{log}" \
        | gzip -c > "{output}"
        """

