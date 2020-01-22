#!/usr/bin/env python

localrules:
    map_counts_to_annotations,
    combine_annotation_counts,
    diffbind_results_to_narrowpeak

rule map_counts_to_annotations:
    input:
        bed = lambda wc: "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_{species}-{factor}-peaks.bed" if wc.annotation=="peaks" else config["differential_occupancy"]["annotations"][wc.annotation],
        bg = lambda wc: "coverage/counts/{sample}_{factor}-chipseq-counts-midpoints.bedgraph" if wc.species=="experimental" else "coverage/sicounts/{sample}_{factor}-chipseq-sicounts-midpoints.bedgraph"
    output:
        temp("diff_binding/{annotation}/{condition}-v-{control}/{sample}_{species}-{factor}-chipseq-counts-{annotation}.tsv")
    log:
        "logs/map_counts_to_annotations/map_counts_to_annotations-{condition}-v-{control}-{sample}-{species}-{annotation}-{factor}.log"
    shell: """
        (cut -f1-6 {input.bed} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule combine_annotation_counts:
    input:
        lambda wc: ["diff_binding/{annotation}/{condition}-v-{control}/".format(**wc) + x + f"_{wc.species}-{wc.factor}-chipseq-counts-{wc.annotation}.tsv" for x in get_samples(search_dict=SAMPLES, passing=True, groups=[wc.control, wc.condition])]
    output:
        "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-{species}-{factor}-chipseq-counts-{annotation}.tsv.gz"
    params:
        n = lambda wc: 7*len(get_samples(search_dict=SAMPLES, passing=True, groups=[wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples(search_dict=SAMPLES, passing=True, groups=[wc.control, wc.condition]).keys())
    log:
        "logs/combine_transcript_counts/combine_transcript_counts-{condition}-v-{control}-{species}-{annotation}-{factor}.log"
    shell: """
        (paste {input} | \
         cut -f$(paste -d, <(echo "1-6") <(seq -s, 7 7 {params.n})) | \
         cat <(echo -e "chrom\tstart\tend\tname\tscore\tstrand\t{params.names}") - | \
         pigz -f > {output}) &> {log}
        """

rule differential_binding:
    input:
        exp_counts = "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-experimental-{factor}-chipseq-counts-{annotation}.tsv.gz",
        spike_counts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{factor}-chipseq-counts-peaks.tsv.gz"
    output:
        counts_norm = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-chipseq-counts-sizefactornorm.tsv",
        counts_rlog = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-chipseq-counts-rlogtransform.tsv",
        results_all = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-all.tsv",
        results_up = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-up.tsv",
        results_down = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-down.tsv",
        results_nonsig = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-nonsignificant.tsv",
        qc_plots = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-qcplots.svg",
    params:
        samples = lambda wc: list(get_samples(search_dict=SAMPLES,
                                              passing=True,
                                              spikein=(True if wc.norm=="spikenorm" else False),
                                              groups=[wc.control, wc.condition]).keys()),
        conditions = lambda wc: [v["group"] for k,v in get_samples(search_dict=SAMPLES,
                                                                  passing=True,
                                                                  spikein=(True if wc.norm=="spikenorm" else False),
                                                                  groups=[wc.control, wc.condition]).items()],
        sampletypes = lambda wc: [("input" if k in INPUTS else "ChIP") \
                                    for k in get_samples(search_dict=SAMPLES,
                                                         passing=True,
                                                         spikein=(True if wc.norm=="spikenorm" else False),
                                                         groups=[wc.control, wc.condition]).keys()],
        alpha = config["differential_occupancy"]["fdr"],
        lfc = log2(config["differential_occupancy"]["fold-change-threshold"])
    conda:
        "../envs/diff_exp.yaml"
    script:
        "../scripts/differential_binding_chipseq.R"

rule diffbind_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand(f"coverage/{wc.norm}/{{sample}}_{FACTOR}-chipseq-{wc.norm}-midpoints-input-subtracted_smoothed.bw",
                sample=get_samples(passing=True,
                                   spikein=(True if wc.norm=="spikenorm" else False),
                                   paired=True,
                                   groups=[wc.condition])),
        control_coverage = lambda wc: expand(f"coverage/{wc.norm}/{{sample}}_{FACTOR}-chipseq-{wc.norm}-midpoints-input-subtracted_smoothed.bw",
                sample=get_samples(passing=True,
                                   spikein=(True if wc.norm=="spikenorm" else False),
                                   paired=True,
                                   groups=[wc.control])),
        diffbind_results = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}.tsv",
    output:
        narrowpeak = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}.narrowpeak",
        summit_bed = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}-summits.bed",
    conda:
        "../envs/diffbind_to_narrowpeak.yaml"
    log:
        "logs/diffbind_results_to_narrowpeak/diffbind_results_to_narrowpeak-{condition}-v-{control}_{annotation}-{norm}-{direction}-{factor}.log"
    shell: """
        (python scripts/chipseq_diffbind_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -d {input.diffbind_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """

