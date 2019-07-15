#!/usr/bin/env python

rule map_to_windows:
    input:
        bg = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-midpoints.bedgraph",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        temp("qual_ctrl/scatter_plots/{factor}-chipseq_{sample}-{norm}-window-{windowsize}.bedgraph")
    log:
        "logs/map_to_windows/map_to_windows_{norm}_{sample}-{windowsize}-{factor}.log"
    shell: """
        (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes) -w {wildcards.windowsize} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule join_window_counts:
    input:
        lambda wc: expand(f"qual_ctrl/scatter_plots/{FACTOR}-chipseq_{{sample}}-{wc.norm}-window-{wc.windowsize}.bedgraph", sample=(SAMPLES if wc.norm=="libsizenorm" else SISAMPLES))
    output:
        "qual_ctrl/scatter_plots/{factor}-chipseq_union-bedgraph-{norm}-window-{windowsize}-allsamples.tsv.gz"
    params:
        names = lambda wc: list(SAMPLES.keys()) if wc.norm=="libsizenorm" else list(SISAMPLES.keys())
    log:
        "logs/join_window_counts/join_window_counts-{norm}-{windowsize}-{factor}.log"
    shell: """
        (bedtools unionbedg -i {input} -header -names {params.names} | \
         bash scripts/cleanUnionbedg.sh | \
         pigz -f > {output}) &> {log}
        """

rule plot_scatter_plots:
    input:
        "qual_ctrl/scatter_plots/{factor}-chipseq_union-bedgraph-{norm}-window-{windowsize}-allsamples.tsv.gz"
    output:
        "qual_ctrl/scatter_plots/{condition}-v-{control}/{status}/{condition}-v-{control}_{factor}-chipseq-{norm}-scatterplots-{status}-window-{windowsize}.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = lambda wc: get_samples(wc.status, wc.norm, [wc.condition, wc.control])
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_scatter_plots.R"

