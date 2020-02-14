#!/usr/bin/env python

localrules:

rule quantify_shift_changes:
    input:
        coverage = lambda wc: expand("coverage/libsizenorm/{sample}_{{factor}}-chipseq-libsizenorm-ratio.bw",
                sample=get_samples(passing=True,
                                   paired=True,
                                   groups=[wc.control, wc.condition])),
        annotation = lambda wc: config["shifts"][wc.annotation]
    output:
        "shifts/{annotation}/{condition}-v-{control}/{condition}-v-{control}_{factor}-chipseq-{annotation}-shifts.tsv"
    params:
        groups = lambda wc: [v["group"] for k,v in get_samples(passing=True,
                                                               paired=True,
                                                               groups=[wc.control, wc.condition]).items()]
    conda:
        "../envs/shifts.yaml"

    shell: """
        python scripts/find_shift_changes.py -i {input.coverage} -g {params.groups} -c {wildcards.control} -b {input.annotation} -o {output}
        """
