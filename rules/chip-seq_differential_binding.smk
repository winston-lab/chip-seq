#!/usr/bin/env python

localrules:
    map_counts_to_annotations,
    combine_annotation_counts

rule map_counts_to_annotations:
    input:
        bed = lambda wc: "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_{species}-{factor}-peaks.bed" if wc.annotation=="peaks" else config["differential_occupancy"]["annotations"][wc.annotation],
        bg = lambda wc: "coverage/counts/{sample}_{factor}-chipseq-counts-midpoints.bedgraph" if wc.species=="experimental" else "coverage/sicounts/{sample}_{factor}-chipseq-sicounts-midpoints.bedgraph"
    output:
        temp("diff_binding/{annotation}/{condition}-v-{control}/{sample}_{species}-{factor}-chipseq-counts-{annotation}.tsv")
    log:
        "logs/map_counts_to_annotations/map_counts_to_annotations-{condition}-v-{control}-{sample}-{species}-{annotation}-{factor}.log"
    shell: """
        (LC_COLLATE=C sort -k1,1 -k2,2n {input.bed} | \
         bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule combine_annotation_counts:
    input:
        lambda wc: ["diff_binding/{annotation}/{condition}-v-{control}/".format(**wc) + x + f"_{wc.species}-{wc.factor}-chipseq-counts-{wc.annotation}.tsv" for x in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-{species}-{factor}-chipseq-counts-{annotation}.tsv.gz"
    params:
        n = lambda wc: 7*len(get_samples("passing", "libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log:
        "logs/combine_transcript_counts/combine_transcript_counts-{condition}-v-{control}-{species}-{annotation}-{factor}.log"
    shell: """
        (paste {input} | \
         cut -f$(paste -d, <(echo "1-6") <(seq -s, 7 7 {params.n})) | \
         cat <(echo -e "chrom\tstart\tend\tname\tscore\tstrand\t{params.names}" ) - > {output}) &> {log}
        """

rule differential_binding:
    input:
        exp_counts = "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-experimental-{factor}-chipseq-counts-{annotation}.tsv.gz",
        spike_counts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{factor}-chipseq-counts-peaks.tsv.gz"
    output:
        counts_norm_chip = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-ChIP-counts-sizefactornorm.tsv",
        counts_norm_input = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-input-counts-sizefactornorm.tsv",
        counts_rlog_chip = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-ChIP-counts-rlogtransform.tsv",
        counts_rlog_input = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-input-counts-rlogtransform.tsv",
        results_all = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-all.tsv",
        results_up = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-up.tsv",
        results_down = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-down.tsv",
        results_unchanged = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-unchanged.tsv",
        bed_all = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-all.bed",
        bed_up = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-up.bed",
        bed_down = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-down.bed",
        bed_unchanged = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-unchanged.bed",
        qc_plots = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-qcplots.svg",
    params:
        chip_samples = lambda wc: [x for x in get_samples("passing", wc.norm, [wc.control, wc.condition]) if x in CHIPS],
        input_samples = lambda wc: [x for x in get_samples("passing", wc.norm, [wc.control, wc.condition]) if x in INPUTS],
        chip_groups = lambda wc: [CHIPS[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition]) if x in CHIPS],
        input_groups = lambda wc: [INPUTS[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition]) if x in INPUTS],
        alpha = config["differential_occupancy"]["fdr"],
        lfc = log2(config["differential_occupancy"]["fold-change-threshold"])
    conda:
        "../envs/diff_exp.yaml"
    script:
        "../scripts/differential_binding.R"

