#!/usr/bin/env python

localrules:
    map_counts_to_annotations,
    combine_annotation_counts,

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
        "../envs/diff_binding.yaml"
    script:
        "../scripts/differential_binding_chipseq_factornorm.R"


