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

include: "rules/differential_binding.smk"
include: "rules/genome_coverage.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: target

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
        expand("coverage/ratio_coverage/libsizenorm/{group}_{nfactor}-over-{dfactor}_libsizenorm-ratio-coverage-window-{windowsize}.bw",
                nfactor=FACTORS["numerator"],
                dfactor=FACTORS["denominator"],
                group=set(conditiongroups + controlgroups),
                windowsize=config["coverage_binsize"]) if comparisons else [],
        expand("coverage/ratio_coverage/spikenorm/{group}_{nfactor}-over-{dfactor}_spikenorm-ratio-coverage-window-{windowsize}.bw",
                nfactor=FACTORS["numerator"],
                dfactor=FACTORS["denominator"],
                group=set(conditiongroups + controlgroups),
                windowsize=config["coverage_binsize"]) if comparisons_si else []

