#!/usr/bin/env python

localrules:
    normalize_genome_coverage,
    # subtract_inputs,
    bedgraph_to_bigwig,
    smoothed_midpoint_coverage,
    map_counts_to_windows,
    combine_window_counts

rule crosslink_coverage:
    input:
        lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-".format(**wc) + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
    output:
        "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-{readtype}.bedgraph",
    params:
        strand_symbol = lambda wc: {"plus": "+", "minus": "-"}.get(wc.readtype)
    wildcard_constraints:
        counttype="counts|sicounts",
        readtype="plus|minus"
    log:
        "logs/crosslink_coverage/crosslink_coverage-{sample}-{counttype}-{readtype}-{factor}.log"
    shell: """
        (bedtools genomecov -bga -5 -strand {params.strand_symbol} -ibam {input} | \
         LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

# extend reads to the median fragment size over all passing samples as
# called by MACS2 cross-correlation
rule protection_coverage:
    input:
        tsv = lambda wc: expand(f"peakcalling/sample_peaks/{{sample}}_{{species}}-{FACTOR}-chipseq_peaks.xls",
                                sample=get_samples(passing=True, paired=True),
                                species=("experimental" if wc.counttype=="counts" else "spikein")),
        bam = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-".format(**wc) + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
    output:
        "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-protection.bedgraph",
    wildcard_constraints:
        counttype="counts|sicounts"
    log:
        "logs/genome_coverage/genome_coverage-{sample}-{counttype}-protection-{factor}.log"
    shell: """
        median_fragsize=$(grep -e "^# d = " {input.tsv} | \
                          cut -d ' ' -f4 | \
                          sort -k1,1n | \
                          awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 2.0;}} }}' | \
                          xargs printf "%.*f\n" 0)
        (bedtools genomecov -bga -fs $median_fragsize -scale $(echo 1/$median_fragsize | bc -l) -ibam {input.bam} | \
         LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

#shift 5' ends of reads by half the median fragment size
#over all samples
rule midpoint_coverage:
    input:
        tsv = lambda wc: expand(f"peakcalling/sample_peaks/{{sample}}_{{species}}-{FACTOR}-chipseq_peaks.xls",
                                sample=get_samples(passing=True, paired=True),
                                species=("experimental" if wc.counttype=="counts" else "spikein")),
        fasta = lambda wc: os.path.abspath(build_annotations(config["genome"]["fasta"])) if wc.counttype=="counts" else config["spike_in"]["fasta"],
        plus = "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-minus.bedgraph"
    output:
        "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-midpoints.bedgraph",
    wildcard_constraints:
        counttype="counts|sicounts"
    log:
        "logs/genome_coverage/genome_coverage-{sample}-{counttype}-midpoints-{factor}.log"
    shell: """
        half_median_fragsize=$(grep -e "^# d = " {input.tsv} | \
                               cut -d ' ' -f4 | \
                               sort -k1,1n | \
                               awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]/2.0}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 4.0;}} }}' | \
                               xargs printf "%.*f\n" 0)
        (bedtools unionbedg -i <(bedtools shift -i {input.plus} -g <(faidx {input.fasta} -i chromsizes) -s $half_median_fragsize | bedtools groupby -g 1,2,3 -c 4 -o sum) <(bedtools shift -i {input.minus} -g <(faidx {input.fasta} -i chromsizes) -s -$half_median_fragsize | bedtools groupby -g 1,2,3 -c 4 -o sum) -g <(faidx {input.fasta} -i chromsizes) -empty | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4+$5}}' > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_{factor}-chipseq-counts-{readtype}.bedgraph",
        bam_experimental = "alignment/{sample}_{factor}-chipseq-uniquemappers-experimental.bam",
        bam_spikein = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-spikein.bam" if wc.norm=="spikenorm" and wc.sample in CHIPS else [],
        input_bam_experimental = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-experimental.bam".format(sample=CHIPS[wc.sample]["control"], factor=wc.factor) if wc.norm=="spikenorm" and wc.sample in CHIPS else [],
        input_bam_spikein = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-spikein.bam".format(sample=CHIPS[wc.sample]["control"], factor=wc.factor) if wc.norm=="spikenorm" and wc.sample in CHIPS else []
    output:
        normalized = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{readtype}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        readtype="plus|minus|protection|midpoints"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{readtype}-{factor}.log"
    run:
        if wildcards.norm=="libsizenorm" or wildcards.sample in INPUTS:
            shell("""
                  (awk -v norm_factor=$(samtools view -c {input.bam_experimental} | \
                                        paste -d "" - <(echo "/1000000") | bc -l) \
                        'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
                  """)
        else:
            shell("""
                  (awk -v norm_factor=$(paste -d "" \
                          <(samtools view -c {input.bam_spikein}) <(echo "*") \
                          <(samtools view -c {input.input_bam_experimental}) <(echo "/") \
                          <(samtools view -c {input.input_bam_spikein}) <(echo "/1000000") | bc -l) \
                          'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
                  """)

rule map_counts_to_windows:
    input:
        bedgraph = "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-midpoints.bedgraph",
        fasta = lambda wc: {"counts": os.path.abspath(build_annotations(config["genome"]["fasta"])),
                            "sicounts": config["spike_in"]["fasta"]
                            }.get(wc.counttype)
    output:
        temp("coverage/{counttype}/{factor}_chipseq_{sample}-{counttype}-midpoints-window-{windowsize}.bedgraph")
    log:
        "logs/map_to_windows/map_to_windows_{sample}-{factor}-{counttype}-{windowsize}.log"
    shell: """
        (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes) -w {wildcards.windowsize} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         bedtools map -a stdin -b {input.bedgraph} -c 4 -o sum > {output}) &> {log}
        """

rule combine_window_counts:
    input:
        bedgraphs = lambda wc: expand("coverage/{{counttype}}/{{factor}}_chipseq_{sample}-{{counttype}}-midpoints-window-{{windowsize}}.bedgraph",
                          sample=get_samples(search_dict=SAMPLES,
                                             passing=True,
                                             groups=[wc.group])),
        fasta = lambda wc: {"counts": os.path.abspath(build_annotations(config["genome"]["fasta"])),
                            "sicounts": config["spike_in"]["fasta"]
                            }.get(wc.counttype)
    output:
        "coverage/ratio_coverage/{factor}_chipseq_{group}-{counttype}-midpoints-window-{windowsize}.tsv.gz"
    params:
        names = lambda wc: list(get_samples(search_dict=SAMPLES,
                                            passing=True,
                                            groups=[wc.group]).keys())
    log:
        "logs/join_window_counts/join_window_counts_{factor}-{group}-{counttype}-{windowsize}.log"
    shell: """
        (bedtools unionbedg -i {input.bedgraphs} -g <(faidx -i chromsizes {input.fasta}) -empty -header -names {params.names} | \
            bash scripts/cleanUnionbedg.sh | \
            pigz -f > {output}) &> {log}
        """

rule ratio_coverage:
    input:
        exp_table = "coverage/ratio_coverage/{factor}_chipseq_{group}-counts-midpoints-window-{windowsize}.tsv.gz",
        spike_table = lambda wc: [] if wc.norm=="libsizenorm" else "coverage/ratio_coverage/{factor}_chipseq_{group}-sicounts-midpoints-window-{windowsize}.tsv.gz",
    output:
        counts_norm = "coverage/ratio_coverage/{norm}/{group}_{factor}_{norm}-ratio-coverage-counts-window-{windowsize}-sizefactornorm.tsv.gz",
        counts_rlog = "coverage/ratio_coverage/{norm}/{group}_{factor}_{norm}-ratio-coverage-counts-window-{windowsize}-rlogtransform.tsv.gz",
        tsv = "coverage/ratio_coverage/{norm}/{group}_{factor}_{norm}-ratio-coverage-counts-window-{windowsize}-results.tsv.gz",
        qc_plots = "coverage/ratio_coverage/{norm}/{group}_{factor}_{norm}-ratio-coverage-counts-window-{windowsize}-qcplots.svg",
        bedgraph = "coverage/{norm}/{group}_{factor}-chipseq-{norm}-ratio_window_{windowsize}.bedgraph",
    params:
        samples = lambda wc: list(get_samples(search_dict=SAMPLES,
                                              passing=True,
                                              spikein=(wc.norm=="spikenorm"),
                                              groups=[wc.group]).keys()),
        rna_sources = lambda wc: [("input" if k in INPUTS else "ChIP") \
                for k in get_samples(search_dict=SAMPLES,
                                     passing=True,
                                     spikein=(wc.norm=="spikenorm"),
                                     groups=[wc.group]).keys()],
    conda:
        "../envs/diff_exp.yaml"
    script:
        "../scripts/chipseq_shrunken_ratio_coverage.R"

rule bedgraph_to_bigwig:
    input:
        bg = "coverage/{norm}/{sample_group}_{factor}-chipseq-{norm}-{readtype}.bedgraph",
        fasta = lambda wc: os.path.abspath(config["spike_in"]["fasta"]) if wc.norm=="sicounts" else os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "coverage/{norm}/{sample_group}_{factor}-chipseq-{norm}-{readtype}.bw"
    log :
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample_group}-{norm}-{readtype}-{factor}.log"
    shell: """
        (bedGraphToBigWig {input.bg} <(faidx {input.fasta} -i chromsizes) {output}) &> {log}
        """

rule smoothed_midpoint_coverage:
    input:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bw"
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}_smoothed.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    conda:
        "../envs/smooth_coverage.yaml"
    log:
        "logs/smoothed_midpoint_coverage/smoothed_midpoint_coverage_{sample}-{norm}-{readtype}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """

