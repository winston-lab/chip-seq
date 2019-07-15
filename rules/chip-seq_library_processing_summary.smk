#!/usr/bin/env python

localrules:
    aggregate_read_numbers,
    plot_read_processing,
    build_spikein_counts_table,
    plot_spikein_pct,

rule aggregate_read_numbers:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("logs/align/align_{sample}.log", sample=SAMPLES),
    output:
        f"qual_ctrl/read_processing/{FACTOR}-chipseq_read-processing-summary.tsv"
    log:
        "logs/aggregate_read_numbers.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map" > {output}) &> {log}""")
        for sample, adapter, align in zip(SAMPLES.keys(), input.adapter, input.align):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | \
                      cut -d: -f2 | \
                      sed 's/,//g' | \
                      awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &>> {log}""")
            shell("""(grep -e "aligned 0 times" -e "aligned exactly 1 time" {align} | \
                      awk 'BEGIN{{ORS="\t"}} {{print $1}} END{{ORS="\\n"; print ""}}' >> {output}) &>> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$4=$3-$4; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &>> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing/{factor}-chipseq_read-processing-summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing/{factor}-chipseq_read-processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing/{factor}-chipseq_read-processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing/{factor}-chipseq_read-processing-loss.svg",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/processing_summary.R"

rule build_spikein_counts_table:
    input:
        ip_bam_experimental = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-experimental.bam", sample=get_samples(spikein=True, paired=True)),
        ip_bam_spikein = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-spikein.bam", sample=get_samples(spikein=True, paired=True)),
        input_bam_experimental = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-experimental.bam", sample=[v["control"] for k,v in get_samples(spikein=True, paired=True).items()]),
        input_bam_spikein = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-spikein.bam", sample=[v["control"] for k,v in get_samples(spikein=True, paired=True).items()])
    output:
        f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-counts.tsv"
    params:
        groups = [v["group"] for k,v in get_samples(spikein=True, paired=True).items()]
    log:
        "logs/build_spikein_counts_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_counts_input\texperimental_counts_input\tspikein_counts_input\ttotal_counts_IP\texperimental_counts_IP\tspikein_counts_IP" > {output}) &> {log} """)
        for sample, group, input_exp, input_si ,ip_exp, ip_si in zip(get_samples(spikein=True, paired=True).keys(), params.groups,
                                                                     input.input_bam_experimental, input.input_bam_spikein,
                                                                     input.ip_bam_experimental, input.ip_bam_spikein):
            shell("""(paste <(echo -e "{sample}\t{group}\t") \
                        <(samtools view -c {input_exp}) \
                        <(samtools view -c {input_si}) \
                        <(echo "") \
                        <(samtools view -c {ip_exp}) \
                        <(samtools view -c {ip_si}) | \
                        awk 'BEGIN{{FS=OFS="\t"}} {{$3=$4+$5; $6=$7+$8; print $0}}'>> {output}) &>> {log} """)

rule plot_spikein_pct:
    input:
        f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-counts.tsv"
    output:
        plot = f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-plots-{{status}}.svg",
        stats = f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-stats-{{status}}.tsv"
    params:
        samplelist = lambda wc: get_samples(passing=(True if wc.status=="passing" else False), spikein=True, paired=True).keys(),
        conditions = conditiongroups_si if comparisons_si else [],
        controls = controlgroups_si if comparisons_si else []
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/spikein_abundance_chipseq.R"

