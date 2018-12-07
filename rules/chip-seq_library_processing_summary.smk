#!/usr/bin/env python

localrules:
    aggregate_read_numbers,
    build_spikein_counts_table,
    plot_spikein_pct,

rule aggregate_read_numbers:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("logs/align/align_{sample}.log", sample=SAMPLES),
    output:
        f"qual_ctrl/read_processing/{FACTOR}-chipseq_read-processing-summary.tsv"
    log: "logs/aggregate_read_numbers.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map" > {output}) &> {log}""")
        for sample, adapter, align in zip(SAMPLES.keys(), input.adapter, input.align):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &>> {log}""")
            shell("""(grep -e "aligned 0 times" -e "aligned exactly 1 time" {align} | awk 'BEGIN{{ORS="\t"}} {{print $1}}' >> {output}) &>> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$4=$3-$4; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &>> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing/{factor}-chipseq_read-processing-summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing/{factor}-chipseq_read-processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing/{factor}-chipseq_read-processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing/{factor}-chipseq_read-processing-loss.svg",
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/processing_summary.R"

rule build_spikein_counts_table:
    input:
        # TODO: handle spikeins for input samples
        # total_bam = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam", sample=SISAMPLES),
        # exp_bam = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-experimental.bam", sample=SISAMPLES),
        # si_bam = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-spikein.bam", sample=SISAMPLES),
        total_bam = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam", sample=CHIPS_SISAMPLES),
        exp_bam = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-experimental.bam", sample=CHIPS_SISAMPLES),
        si_bam = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-spikein.bam", sample=CHIPS_SISAMPLES),
    output:
        f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-counts.tsv"
    params:
        # groups = [v["group"] for k,v in SISAMPLES.items()]
        groups = [v["group"] for k,v in CHIPS_SISAMPLES.items()]
    log: "logs/build_spikein_counts_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_counts\texperimental_counts\tspikein_counts" > {output}) &> {log} """)
        # for sample, group, total, exp, si in zip(SISAMPLES.keys(), params.groups, input.total_bam, input.exp_bam, input.si_bam):
        for sample, group, total, exp, si in zip(CHIPS_SISAMPLES.keys(), params.groups, input.total_bam, input.exp_bam, input.si_bam):
            shell("""(paste <(echo -e "{sample}\t{group}") <(samtools view -c {total}) <(samtools view -c {exp}) <(samtools view -c {si}) >> {output}) &>> {log} """)

rule plot_spikein_pct:
    input:
        f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-counts.tsv"
    output:
        plot = f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-plots-{{status}}.svg",
        stats = f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-stats-{{status}}.tsv"
    params:
        # samplelist = lambda wc : list(SISAMPLES.keys()) if wc.status=="all" else list(SIPASSING.keys()),
        samplelist = lambda wc : list(CHIPS_SISAMPLES.keys()) if wc.status=="all" else list(CHIPS_SIPASSING.keys()),
        conditions = conditiongroups_si if comparisons_si else [],
        controls = controlgroups_si if comparisons_si else []
    conda: "../envs/tidyverse.yaml"
    script: "../scripts/plot_si_pct.R"

