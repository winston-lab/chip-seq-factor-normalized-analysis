#!/usr/bin/env python

localrules:
    quantify_shift_changes

rule quantify_shift_changes:
    input:
        coverage = lambda wc: expand("coverage/ratio_coverage/libsizenorm/{group}_{{nfactor}}-over-{{dfactor}}_libsizenorm-standard-difference.bw",
                group=[wc.control, wc.condition]),
        annotation = lambda wc: config["shifts"][wc.annotation]
    output:
        "shifts/{annotation}/{condition}-v-{control}/{condition}-v-{control}_{nfactor}-over-{dfactor}-chipseq-{annotation}-shifts.tsv"
    conda:
        "../envs/shifts.yaml"
    shell: """
        python scripts/find_shift_changes_factornorm.py -i {input.coverage} -g {wildcards.control} {wildcards.condition} -c {wildcards.control} -b {input.annotation} -o {output}
        """
