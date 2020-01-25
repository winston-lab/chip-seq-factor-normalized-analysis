#!/usr/bin/env python

localrules:
    compute_matrix,
    cat_matrices

rule compute_matrix:
    input:
        annotation = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["path"],
        bw = lambda wc: f"coverage/ratio_coverage/{{norm}}/{{sample}}_{FACTORS['numerator']}-over-{FACTORS['denominator']}_{{norm}}-ratio-coverage-window-{config['coverage_binsize']}.bw" if wc.strand == "ratio" else denominator_pipe(f"coverage/{{norm}}/{{sample}}_{FACTORS['denominator']}-chipseq-{{norm}}-{{strand}}.bw"),
    output:
        dtfile = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}.tsv"),
        melted = temp("datavis/{figure}/{norm}/{annotation}_{sample}-{norm}-{strand}-melted.tsv.gz"),
    params:
        group = lambda wc: wc.sample if wc.strand == "ratio" else SAMPLES["denominator"][wc.sample]["group"],
        refpoint = lambda wc: "TSS" if FIGURES[wc.figure]["parameters"]["type"]=="scaled" else FIGURES[wc.figure]["parameters"]["refpoint"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        binsize = lambda wc: FIGURES[wc.figure]["parameters"]["binsize"],
        binstat = lambda wc: FIGURES[wc.figure]["parameters"]["binstat"],
        nan_afterend = lambda wc: [] if FIGURES[wc.figure]["parameters"]["type"]=="scaled" or not FIGURES[wc.figure]["parameters"]["nan_afterend"] else "--nanAfterEnd",
        anno_label = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["label"]
    threads:
        config["threads"]
    log:
        "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}-{norm}-{strand}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} --group {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{{figure}}/{{norm}}/{annotation}_{sample}-{{norm}}-{{strand}}-melted.tsv.gz",
                          annotation=list(FIGURES[wc.figure]["annotations"].keys()),
                          sample={"libsizenorm": set(conditiongroups + controlgroups),
                                  "spikenorm": set(conditiongroups_si + controlgroups_si)}.get(wc.norm) if wc.strand=="ratio" else \
                                  get_samples(paired=True,
                                              passing=True,
                                              spikein=(wc.norm == "spikenorm")))
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-chipseq-{norm}-{strand}.tsv.gz"
    log:
        "logs/cat_matrices/cat_matrices-{figure}_{norm}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()],
        denominator = "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-chipseq-{norm}-{strand}.tsv.gz",
        ratio = "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-chipseq-{norm}-ratio.tsv.gz"
    output:
        ratio_heatmap = "datavis/{figure}/{norm}/{condition}-v-{control}/{strand}/{figure}_{condition}-v-{control}_{nfactor}-over-{dfactor}_{norm}-ratio-heatmap.svg",
        ratio_metagene = "datavis/{figure}/{norm}/{condition}-v-{control}/{strand}/{figure}_{condition}-v-{control}_{nfactor}-over-{dfactor}_{norm}-ratio-metagene.svg",
        ridgelines = "datavis/{figure}/{norm}/{condition}-v-{control}/{strand}/{figure}_{condition}-v-{control}_{nfactor}-over-{dfactor}_{norm}-ratio-{strand}-ridgelines.svg"
    params:
        # abusing snakemake a bit here...using params as output paths since in order to use lambda functions
        annotations_out = lambda wc: expand(f"datavis/{wc.figure}/{wc.norm}/{wc.condition}-v-{wc.control}/{wc.strand}/{{annotation}}.bed", annotation=FIGURES[wc.figure]["annotations"]),
        length_sort = lambda wc: FIGURES[wc.figure]["parameters"]["length_sort"],
        distance_type = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        colorbar_quantile_low = lambda wc: FIGURES[wc.figure]["parameters"]["colorbar_quantile_low"],
        colorbar_quantile_high = lambda wc: FIGURES[wc.figure]["parameters"]["colorbar_quantile_high"],
        refpt_label = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        endpt_label = lambda wc: "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["endlabel"],
        heatmap_colormap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        ridgeline_colormap = lambda wc: FIGURES[wc.figure]["parameters"]["ridgeline_colormap"],
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/chipseq_factornorm_datavis.R"

