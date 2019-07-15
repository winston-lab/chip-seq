#!/usr/bin/env python

rule map_to_windows:
    input:
        bg = f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-midpoints.bedgraph",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        temp(f"qual_ctrl/scatter_plots/{FACTOR}_chipseq_{{sample}}-{{norm}}-midpoints-window-{{windowsize}}.bedgraph")
    log:
        "logs/map_to_windows/map_to_windows_{sample}-{norm}-{windowsize}.log"
    shell: """
        (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes) -w {wildcards.windowsize} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule join_window_counts:
    input:
        lambda wc: expand("qual_ctrl/scatter_plots/{factor}_chipseq_{sample}-{{norm}}-midpoints-window-{{windowsize}}.bedgraph", sample=(SAMPLES if wc.norm=="libsizenorm" else SISAMPLES), factor=FACTOR)
    output:
        f"qual_ctrl/scatter_plots/{FACTOR}_chipseq_union-bedgraph-{{norm}}-midpoint-window-{{windowsize}}-allsamples.tsv.gz"
    params:
        names = lambda wc: list((SAMPLES if wc.norm=="libsizenorm" else SISAMPLES).keys())
    log:
        "logs/join_window_counts/join_window_counts-{norm}-{windowsize}.log"
    shell: """
        (bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output}) &> {log}
        """

rule plot_scatter_plots:
    input:
        "qual_ctrl/scatter_plots/{factor}_chipseq_union-bedgraph-{norm}-midpoint-window-{windowsize}-allsamples.tsv.gz"
    output:
        "qual_ctrl/scatter_plots/{condition}-v-{control}/{status}/{condition}-v-{control}_{factor}_chipseq-{norm}-scatterplots-{status}-window-{windowsize}.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = lambda wc: list(get_samples(search_dict=SAMPLES,
                                                 passing=(True if wc.status=="passing" else False),
                                                 spikein=(True if wc.norm=="spikenorm" else False),
                                                 groups=[wc.condition, wc.control]).keys())
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_scatter_plots.R"

