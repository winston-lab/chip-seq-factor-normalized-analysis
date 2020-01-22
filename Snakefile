#!/usr/bin/env python

import yaml
import itertools

configfile: "config.yaml"

subworkflow denominator_pipe:
    workdir: config["denominator_pipe"]
subworkflow numerator_pipe:
    workdir: config["numerator_pipe"]

with open(denominator_pipe("config.yaml"), "r") as stream:
    config_d = yaml.safe_load(stream)
with open(numerator_pipe("config.yaml"), "r") as stream:
    config_n = yaml.safe_load(stream)

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

def get_samples(search_dict=CHIPS["denominator"],
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
                if v["control"] in get_samples(search_dict=INPUTS["denominator"],
                                               passing=passing,
                                               spikein=spikein,
                                               paired=False)}
    if groups and "all" not in groups:
        search_dict = {k:v for k,v in search_dict.items() if v["group"] in groups}
    return search_dict

rule target:
    input:
        expand(expand("diff_binding/{{annotation}}/{condition}-v-{control}/{condition}-v-{control}_allsamples-{{species}}-{{factor}}-chipseq-counts-{{annotation}}.tsv.gz", zip, condition=conditiongroups, control=controlgroups),
                species=["experimental", "spikein"],
                factor=list(FACTORS.values()),
                annotation=list(config["differential_occupancy"]["annotations"].keys()))

rule map_counts_to_annotations:
    input:
        bed = lambda wc: config["differential_occupancy"]["annotations"][wc.annotation],
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
        lambda wc: ["diff_binding/{annotation}/{condition}-v-{control}/".format(**wc) + x + "_{species}-{factor}-chipseq-counts-{annotation}.tsv".format(**wc) for x in get_samples(search_dict=SAMPLES["denominator" if wc.factor == FACTORS["denominator"] else "numerator"], passing=True, groups=[wc.control, wc.condition])]
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
