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

# rule map_counts_to_windows:
#     input:
#         bedgraph = lambda wc: {True: denominator_pipe("coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-midpoints.bedgraph"),
#                                False: numerator_pipe("coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-midpoints.bedgraph")}.get(wc.factor == FACTORS["denominator"]),
#         fasta = lambda wc: {"counts": os.path.abspath(annotation_pipe(config_annotation["genome"]["fasta"])),
#                             "sicounts": os.path.abspath(denominator_pipe(config_d["spike_in"]["fasta"]))
#                             }.get(wc.counttype)
#     output:
#         temp("coverage/{counttype}/{factor}_chipseq_{sample}-{counttype}-midpoints-window-{windowsize}.bedgraph")
#     log:
#         "logs/map_to_windows/map_to_windows_{sample}-{factor}-{counttype}-{windowsize}.log"
#     shell: """
#         (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes) -w {wildcards.windowsize} | \
#          LC_COLLATE=C sort -k1,1 -k2,2n | \
#          bedtools map -a stdin -b {input.bedgraph} -c 4 -o sum > {output}) &> {log}
#         """

# rule combine_window_counts:
#     input:
#         bedgraphs = lambda wc: expand("coverage/{{counttype}}/{{factor}}_chipseq_{sample}-{{counttype}}-midpoints-window-{{windowsize}}.bedgraph",
#                           sample=get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"],
#                                              passing=True)),
#         fasta = lambda wc: {"counts": os.path.abspath(annotation_pipe(config_annotation["genome"]["fasta"])),
#                             "sicounts": os.path.abspath(denominator_pipe(config_d["spike_in"]["fasta"]))
#                             }.get(wc.counttype)
#     output:
#         "coverage/{counttype}/{factor}_chipseq_allsamples-{counttype}-midpoints-window-{windowsize}.tsv.gz"
#     params:
#         names = lambda wc: list(get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"],
#                                             passing=True).keys())
#     log:
#         "logs/join_window_counts/join_window_counts_{factor}-{counttype}-{windowsize}.log"
#     shell: """
#         (bedtools unionbedg -i {input.bedgraphs} -g <(faidx -i chromsizes {input.fasta}) -empty -header -names {params.names} | \
#             bash scripts/cleanUnionbedg.sh | \
#             pigz -f > {output}) &> {log}
#         """

# rule ratio_coverage:
#     input:
#         exp_table_denominator = "coverage/counts/{dfactor}_chipseq_allsamples-counts-midpoints-window-{windowsize}.tsv.gz",
#         exp_table_numerator = "coverage/counts/{nfactor}_chipseq_allsamples-counts-midpoints-window-{windowsize}.tsv.gz",
#         spike_table_denominator = lambda wc: [] if wc.norm=="libsizenorm" else "coverage/sicounts/{dfactor}_chipseq_allsamples-sicounts-midpoints-window-{windowsize}.tsv.gz",
#         spike_table_numerator = lambda wc: [] if wc.norm=="libsizenorm" else "coverage/sicounts/{nfactor}_chipseq_allsamples-sicounts-midpoints-window-{windowsize}.tsv.gz",
#     output:
#         counts_norm = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-counts-window-{windowsize}-sizefactornorm.tsv.gz",
#         counts_rlog = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-counts-window-{windowsize}-rlogtransform.tsv.gz",
#         tsv = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}-results.tsv.gz",
#         bedgraph = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}.bedgraph",
#         qc_plots = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}-qcplots.svg",
#     params:
#         samples_denominator = lambda wc: list(get_samples(search_dict=SAMPLES["denominator"],
#                                                           passing=True,
#                                                           spikein=(wc.norm=="spikenorm"),
#                                                           groups=[wc.group]).keys()),
#         samples_numerator = lambda wc: list(get_samples(search_dict=SAMPLES["numerator"],
#                                                         passing=True,
#                                                         spikein=(wc.norm=="spikenorm"),
#                                                         groups=[wc.group]).keys()),
#         rna_sources_denominator = lambda wc: [("input" if k in INPUTS["denominator"] else "ChIP") \
#                 for k in get_samples(search_dict=SAMPLES["denominator"],
#                                      passing=True,
#                                      spikein=(wc.norm=="spikenorm"),
#                                      groups=[wc.group]).keys()],
#         rna_sources_numerator = lambda wc: [("input" if k in INPUTS["numerator"] else "ChIP") \
#                 for k in get_samples(search_dict=SAMPLES["numerator"],
#                                      passing=True,
#                                      spikein=(wc.norm=="spikenorm"),
#                                      groups=[wc.group]).keys()],
#     conda:
#         "../envs/diff_binding.yaml"
#     script:
#         "../scripts/chipseq_shrunken_ratio_coverage.R"

# rule bedgraph_to_bigwig:
#     input:
#         bedgraph = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}.bedgraph",
#         fasta = lambda wc: os.path.abspath(annotation_pipe(config_annotation["genome"]["fasta"]))
#     output:
#         "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}.bw",
#     log :
#         "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{group}_{nfactor}-over-{dfactor}-{norm}-{windowsize}.log"
#     shell: """
#         (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes) {output}) &> {log}
#         """

