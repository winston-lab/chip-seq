#!/usr/bin/env python

localrules:
    map_counts_to_annotations,
    combine_transcript_counts

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

rule differential_expression:
    input:
        expcounts = "diff_exp/{condition}-v-{control}/{condition}-v-{control}_{assay}-experimental-transcript-counts.tsv",
        sicounts = lambda wc: [] if wc.norm=="libsizenorm" else  "diff_exp/{condition}-v-{control}/{condition}-v-{control}_{assay}-spikein-transcript-counts.tsv".format(**wc)
    output:
        results_all = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-all.tsv",
        results_up = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-up.tsv",
        results_down = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-down.tsv",
        results_unch = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-unchanged.tsv",
        bed_all = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-all.bed",
        bed_up = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-up.bed",
        bed_down = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-down.bed",
        bed_unch = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-unchanged.bed",
        normcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-counts-sizefactornorm.tsv",
        rldcounts = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-counts-rlogtransformed.tsv",
        qcplots = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-qc-plots.svg"
    params:
        samples = lambda wc : get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc : [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["differential_occupancy"]["fdr"],
        lfc = log2(config["differential_occupancy"]["fold-change-threshold"])
    conda: "../envs/diff_exp.yaml"
    script:
        "../scripts/call_diffexp_transcripts.R"

rule summarise_diffexp_results:
    input:
        total = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-results-all.tsv",
        genic = "diff_exp/{condition}-v-{control}/{norm}/genic/{condition}-v-{control}_{assay}-{norm}-diffexp-results-genic-all.tsv",
        antisense = "diff_exp/{condition}-v-{control}/{norm}/antisense/{condition}-v-{control}_{assay}-{norm}-diffexp-results-antisense-all.tsv",
        convergent = "diff_exp/{condition}-v-{control}/{norm}/convergent/{condition}-v-{control}_{assay}-{norm}-diffexp-results-convergent-all.tsv",
        divergent = "diff_exp/{condition}-v-{control}/{norm}/divergent/{condition}-v-{control}_{assay}-{norm}-diffexp-results-divergent-all.tsv",
        intergenic = "diff_exp/{condition}-v-{control}/{norm}/intergenic/{condition}-v-{control}_{assay}-{norm}-diffexp-results-intergenic-all.tsv",
    output:
        summary_table = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-summary.tsv",
        mosaic = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-mosaic.svg",
        maplot = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-maplot.svg",
        volcano = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-volcano.svg",
        volcano_free = "diff_exp/{condition}-v-{control}/{norm}/{condition}-v-{control}_{assay}-{norm}-diffexp-volcano-freescale.svg",
    params:
        lfc = config["differential_occupancy"]["fold-change-threshold"],
        alpha = config["differential_occupancy"]["fdr"]
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/plot_diffexp_summary.R"

