#!/usr/bin/env python

import os
import yaml
import itertools
from math import log2

configfile: "config.yaml"

subworkflow denominator_pipe:
    workdir: config["denominator_pipe"]
subworkflow numerator_pipe:
    workdir: config["numerator_pipe"]

with open(denominator_pipe("config.yaml"), "r") as stream:
    config_d = yaml.safe_load(stream)
with open(numerator_pipe("config.yaml"), "r") as stream:
    config_n = yaml.safe_load(stream)

annotation_workflow_d = os.path.abspath(config_d["genome"]["annotation_workflow"])
annotation_workflow_n = os.path.abspath(config_n["genome"]["annotation_workflow"])

assert annotation_workflow_d == annotation_workflow_n, "Numerator and denominator annotation workflows do not match."

subworkflow annotation_pipe:
    workdir: annotation_workflow_d

with open(annotation_pipe("config.yaml"), "r") as stream:
    config_annotation = yaml.safe_load(stream)

FACTORS = { "denominator": config_d["factor"],
            "numerator": config_n["factor"]}
assert FACTORS["denominator"] != FACTORS["numerator"], "Numerator factor and denominator factor cannot have the same name."

INPUTS = {  "denominator": config_d["input_samples"],
            "numerator": config_n["input_samples"]}

CHIPS = {   "denominator": config_d["chip_samples"],
            "numerator": config_n["chip_samples"]}

SAMPLES = { "denominator": {**INPUTS["denominator"], **CHIPS["denominator"]},
            "numerator": {**INPUTS["numerator"], **CHIPS["numerator"]}}

comparisons = config["comparisons"]["libsizenorm"]
if comparisons:
    controlgroups = list(itertools.chain(*[d.values() for d in comparisons]))
    conditiongroups = list(itertools.chain(*[d.keys() for d in comparisons]))
comparisons_si = config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in comparisons_si]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in comparisons_si]))
    assert config_d["spike_in"]["fasta"] == config_n["spike_in"]["fasta"], "Spike-in normalized comparisons are specified, but numerator and denominator spike-in FASTA files do not match."

def get_samples(search_dict=CHIPS["denominator"],
                paired_search_dict=None,
                passing=False,
                spikein=False,
                paired=False,
                groups=None):
    if passing:
        search_dict = {k:v for k,v in search_dict.items() if v["pass-qc"]}
    if spikein:
        search_dict = {k:v for k,v in search_dict.items() if v["spikein"]}
    if paired:
        search_dict = {k:v for k,v in search_dict.items() \
                if v["control"] in get_samples(search_dict=paired_search_dict,
                                               passing=passing,
                                               spikein=spikein,
                                               paired=False)}
    if groups and "all" not in groups:
        search_dict = {k:v for k,v in search_dict.items() if v["group"] in groups}
    return search_dict

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules:
    target,
    map_counts_to_annotations,
    combine_annotation_counts,
    map_counts_to_windows,
    combine_window_counts,

rule target:
    input:
        expand(expand("diff_binding/{{annotation}}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{nfactor}}-over-{{dfactor}}-chipseq-libsizenorm-{{annotation}}-diffbind-results-all.tsv",
            zip, condition=conditiongroups, control=controlgroups),
            annotation=list(config["differential_occupancy"]["annotations"].keys()),
            nfactor=FACTORS["numerator"],
            dfactor=FACTORS["denominator"]) if comparisons else [],
        expand(expand("diff_binding/{{annotation}}/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{nfactor}}-over-{{dfactor}}-chipseq-spikenorm-{{annotation}}-diffbind-results-all.tsv",
            zip, condition=conditiongroups_si, control=controlgroups_si),
            annotation=list(config["differential_occupancy"]["annotations"].keys()),
            nfactor=FACTORS["numerator"],
            dfactor=FACTORS["denominator"]) if comparisons_si else [],
        expand("coverage/ratio_coverage/libsizenorm/{group}_{nfactor}-over-{dfactor}_libsizenorm-ratio-coverage-window-{windowsize}.bedgraph",
                nfactor=FACTORS["numerator"],
                dfactor=FACTORS["denominator"],
                group=set(conditiongroups + controlgroups),
                windowsize=config["coverage_binsize"]) if comparisons else [],
        expand("coverage/ratio_coverage/spikenorm/{group}_{nfactor}-over-{dfactor}_spikenorm-ratio-coverage-window-{windowsize}.bedgraph",
                nfactor=FACTORS["numerator"],
                dfactor=FACTORS["denominator"],
                group=set(conditiongroups + controlgroups),
                windowsize=config["coverage_binsize"]) if comparisons_si else []

rule map_counts_to_annotations:
    input:
        bed = lambda wc: config["differential_occupancy"]["annotations"][wc.annotation]["experimental-annotation" if wc.species == "experimental" else "spikein-annotation"],
        bedgraph =
            lambda wc:  {"experimental":
                            {True: denominator_pipe("coverage/counts/{sample}_{factor}-chipseq-counts-midpoints.bedgraph"),
                             False: numerator_pipe("coverage/counts/{sample}_{factor}-chipseq-counts-midpoints.bedgraph")
                            },
                         "spikein":
                            {True: denominator_pipe("coverage/sicounts/{sample}_{factor}-chipseq-sicounts-midpoints.bedgraph"),
                             False: numerator_pipe("coverage/sicounts/{sample}_{factor}-chipseq-sicounts-midpoints.bedgraph")
                            }
                        }.get(wc.species).get(wc.factor == FACTORS["denominator"])
    output:
        temp("diff_binding/{annotation}/{condition}-v-{control}/{sample}_{species}-{factor}-chipseq-counts-{annotation}.tsv")
    log:
        "logs/map_counts_to_annotations/map_counts_to_annotations-{condition}-v-{control}-{sample}-{species}-{annotation}-{factor}.log"
    shell: """
        (cut -f1-6 {input.bed} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bedgraph} -c 4 -o sum > {output}) &> {log}
        """

rule combine_annotation_counts:
    input:
        lambda wc: expand("diff_binding/{{annotation}}/{{condition}}-v-{{control}}/{sample}_{{species}}-{{factor}}-chipseq-counts-{{annotation}}.tsv",
                sample=get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"], passing=True, groups=[wc.control, wc.condition]))
    output:
        "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-{species}-{factor}-chipseq-counts-{annotation}.tsv.gz"
    params:
        n = lambda wc: 7*len(get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"], passing=True, groups=[wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"], passing=True, groups=[wc.control, wc.condition]).keys())
    log:
        "logs/combine_transcript_counts/combine_transcript_counts-{condition}-v-{control}-{species}-{annotation}-{factor}.log"
    shell: """
        (paste {input} | \
         cut -f$(paste -d, <(echo "1-6") <(seq -s, 7 7 {params.n})) | \
         cat <(echo -e "chrom\tstart\tend\tname\tscore\tstrand\t{params.names}" ) - | \
         pigz -f > {output}) &> {log}
        """

rule differential_binding:
    input:
        exp_table_denominator = "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-experimental-{dfactor}-chipseq-counts-{annotation}.tsv.gz",
        exp_table_numerator = "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-experimental-{nfactor}-chipseq-counts-{annotation}.tsv.gz",
        spike_table_denominator = lambda wc: [] if wc.norm=="libsizenorm" else "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{dfactor}-chipseq-counts-{annotation}.tsv.gz",
        spike_table_numerator = lambda wc: [] if wc.norm=="libsizenorm" else "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{nfactor}-chipseq-counts-{annotation}.tsv.gz",
    output:
        counts_norm = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-chipseq-counts-sizefactornorm.tsv",
        counts_rlog = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-chipseq-counts-rlogtransform.tsv",
        results_all = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-all.tsv",
        results_up = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-up.tsv",
        results_down = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-down.tsv",
        results_nonsig = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-nonsignificant.tsv",
        bed_all = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-all.bed",
        bed_up = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-up.bed",
        bed_down = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-down.bed",
        bed_nonsig = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-results-nonsignificant.bed",
        qc_plots = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{norm}-{annotation}-diffbind-qcplots.svg",
    params:
        samples_denominator = lambda wc: list(get_samples(search_dict=SAMPLES["denominator"],
                                                          passing=True,
                                                          spikein=(wc.norm=="spikenorm"),
                                                          groups=[wc.control, wc.condition]).keys()),
        samples_numerator = lambda wc: list(get_samples(search_dict=SAMPLES["numerator"],
                                                        passing=True,
                                                        spikein=(wc.norm=="spikenorm"),
                                                        groups=[wc.control, wc.condition]).keys()),
        conditions_denominator = lambda wc: [v["group"] for k,v in get_samples(search_dict=SAMPLES["denominator"],
                                                                               passing=True,
                                                                               spikein=(wc.norm=="spikenorm"),
                                                                               groups=[wc.control, wc.condition]).items()],
        conditions_numerator = lambda wc: [v["group"] for k,v in get_samples(search_dict=SAMPLES["numerator"],
                                                                             passing=True,
                                                                             spikein=(wc.norm=="spikenorm"),
                                                                             groups=[wc.control, wc.condition]).items()],
        rna_sources_denominator = lambda wc: [("input" if k in INPUTS["denominator"] else "ChIP") \
                for k in get_samples(search_dict=SAMPLES["denominator"],
                                     passing=True,
                                     spikein=(wc.norm=="spikenorm"),
                                     groups=[wc.control, wc.condition]).keys()],
        rna_sources_numerator = lambda wc: [("input" if k in INPUTS["numerator"] else "ChIP") \
                for k in get_samples(search_dict=SAMPLES["numerator"],
                                     passing=True,
                                     spikein=(wc.norm=="spikenorm"),
                                     groups=[wc.control, wc.condition]).keys()],
        alpha = config["differential_occupancy"]["fdr"],
        lfc = log2(config["differential_occupancy"]["fold-change-threshold"])
    conda:
        "envs/diff_binding.yaml"
    script:
        "scripts/differential_binding_chipseq_factornorm.R"

rule map_counts_to_windows:
    input:
        bedgraph = lambda wc: {True: denominator_pipe("coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-midpoints.bedgraph"),
                               False: numerator_pipe("coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-midpoints.bedgraph")}.get(wc.factor == FACTORS["denominator"]),
        fasta = lambda wc: {"counts": os.path.abspath(annotation_pipe(config_annotation["genome"]["fasta"])),
                            "sicounts": os.path.abspath(denominator_pipe(config_d["spike_in"]["fasta"]))
                            }.get(wc.counttype)
    output:
        temp("coverage/{counttype}/{factor}_chipseq_{sample}-{counttype}-midpoints-window-{windowsize}.bedgraph")
    log:
        "logs/map_to_windows/map_to_windows_{sample}-{factor}-{counttype}-{windowsize}.log"
    shell: """
        (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes) -w {wildcards.windowsize} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bedgraph} -c 4 -o sum > {output}) &> {log}
        """

rule combine_window_counts:
    input:
        bedgraphs = lambda wc: expand("coverage/{{counttype}}/{{factor}}_chipseq_{sample}-{{counttype}}-midpoints-window-{{windowsize}}.bedgraph",
                          sample=get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"],
                                             passing=True)),
        fasta = lambda wc: {"counts": os.path.abspath(annotation_pipe(config_annotation["genome"]["fasta"])),
                            "sicounts": os.path.abspath(denominator_pipe(config_d["spike_in"]["fasta"]))
                            }.get(wc.counttype)
    output:
        "coverage/{counttype}/{factor}_chipseq_allsamples-{counttype}-midpoints-window-{windowsize}.tsv.gz"
    params:
        names = lambda wc: list(get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"],
                                            passing=True).keys())
    log:
        "logs/join_window_counts/join_window_counts_{factor}-{counttype}-{windowsize}.log"
    shell: """
        (bedtools unionbedg -i {input.bedgraphs} -g <(faidx -i chromsizes {input.fasta}) -empty -header -names {params.names} | \
            bash scripts/cleanUnionbedg.sh | \
            pigz -f > {output}) &> {log}
        """

rule ratio_coverage:
    input:
        exp_table_denominator = "coverage/counts/{dfactor}_chipseq_allsamples-counts-midpoints-window-{windowsize}.tsv.gz",
        exp_table_numerator = "coverage/counts/{nfactor}_chipseq_allsamples-counts-midpoints-window-{windowsize}.tsv.gz",
        spike_table_denominator = lambda wc: [] if wc.norm=="libsizenorm" else "coverage/sicounts/{dfactor}_chipseq_allsamples-sicounts-midpoints-window-{windowsize}.tsv.gz",
        spike_table_numerator = lambda wc: [] if wc.norm=="libsizenorm" else "coverage/sicounts/{nfactor}_chipseq_allsamples-sicounts-midpoints-window-{windowsize}.tsv.gz",
    output:
        counts_norm = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-counts-window-{windowsize}-sizefactornorm.tsv.gz",
        counts_rlog = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-counts-window-{windowsize}-rlogtransform.tsv.gz",
        tsv = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}-results.tsv.gz",
        bedgraph = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}.bedgraph",
        qc_plots = "coverage/ratio_coverage/{norm}/{group}_{nfactor}-over-{dfactor}_{norm}-ratio-coverage-window-{windowsize}-qcplots.svg",
    params:
        samples_denominator = lambda wc: list(get_samples(search_dict=SAMPLES["denominator"],
                                                          passing=True,
                                                          spikein=(wc.norm=="spikenorm"),
                                                          groups=[wc.group]).keys()),
        samples_numerator = lambda wc: list(get_samples(search_dict=SAMPLES["numerator"],
                                                        passing=True,
                                                        spikein=(wc.norm=="spikenorm"),
                                                        groups=[wc.group]).keys()),
        rna_sources_denominator = lambda wc: [("input" if k in INPUTS["denominator"] else "ChIP") \
                for k in get_samples(search_dict=SAMPLES["denominator"],
                                     passing=True,
                                     spikein=(wc.norm=="spikenorm"),
                                     groups=[wc.group]).keys()],
        rna_sources_numerator = lambda wc: [("input" if k in INPUTS["numerator"] else "ChIP") \
                for k in get_samples(search_dict=SAMPLES["numerator"],
                                     passing=True,
                                     spikein=(wc.norm=="spikenorm"),
                                     groups=[wc.group]).keys()],
    conda:
        "envs/diff_binding.yaml"
    script:
        "scripts/chipseq_shrunken_ratio_coverage.R"

