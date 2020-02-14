#!/usr/bin/env python

localrules:
    factor_ratio_coverage

rule factor_ratio_coverage:
    input:
        numerator_paths = lambda wc: expand(numerator_pipe(f"coverage/{wc.norm}/{{sample}}_{{factor}}-chipseq-{wc.norm}-ratio.bw"),
                sample=get_samples(
                    search_dict=CHIPS["numerator"],
                    paired_search_dict=INPUTS["numerator"],
                    passing=True,
                    spikein=(wc.norm=="spikenorm")),
                factor=FACTORS["numerator"]),
        denominator_paths = lambda wc: expand(denominator_pipe(f"coverage/{wc.norm}/{{sample}}_{{factor}}-chipseq-{wc.norm}-ratio.bw"),
                sample=get_samples(
                    search_dict=CHIPS["denominator"],
                    paired_search_dict=INPUTS["denominator"],
                    passing=True,
                    spikein=(wc.norm=="spikenorm")),
                factor=FACTORS["denominator"])
    output:
        bigwig = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-standard-difference.bw",
    params:
        numerator_groups = lambda wc: [v["group"] for k,v in get_samples(
            search_dict=CHIPS["numerator"],
            paired_search_dict=INPUTS["numerator"],
            passing=True,
            spikein=(wc.norm=="spikenorm")).items()],
        denominator_groups = lambda wc: [v["group"] for k,v in get_samples(
            search_dict=CHIPS["denominator"],
            paired_search_dict=INPUTS["denominator"],
            passing=True,
            spikein=(wc.norm=="spikenorm")).items()]
    conda:
        "../envs/factor_ratio.yaml"
    shell: """
        python scripts/factor_ratio.py -n {input.numerator_paths} -d {input.denominator_paths} -j {params.numerator_groups} -k {params.denominator_groups} -g {wildcards.group} -o {output.bigwig}
        """

