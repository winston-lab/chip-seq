#!/usr/bin/env python

rule crosslink_coverage:
    input:
        lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-".format(**wc) + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
    output:
        "coverage/{counttype}/{sample}_{factor}-chipseq-{counttype}-{strand}.bedgraph",
    params:
        strand_symbol = lambda wc: {"plus": "+", "minus": "-"}.get(wc.strand)
    wildcard_constraints:
        counttype="counts|sicounts",
        strand="plus|minus"
    log:
        "logs/crosslink_coverage/crosslink_coverage-{sample}-{counttype}-{strand}-{factor}.log"
    shell: """
        (bedtools genomecov -bga -5 -strand {params.strand_symbol} -ibam {input} | \
         LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

#extend reads to the median fragment size over all samples as
#called by MACS2 cross-correlation
rule protection_coverage:
    input:
        tsv = lambda wc: expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_peaks.xls", factor=FACTOR, group=GROUPS, species= ("experimental" if wc.counttype=="counts" else "spikein")),
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
        tsv = lambda wc: expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_peaks.xls", factor=FACTOR, group=GROUPS, species= ("experimental" if wc.counttype=="counts" else "spikein")),
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
        (bedtools unionbedg -i <(bedtools shift -i {input.plus} -g <(faidx {input.fasta} -i chromsizes) -s $half_median_fragsize) <(bedtools shift -i {input.minus} -g <(faidx {input.fasta} -i chromsizes) -s -$half_median_fragsize) -g <(faidx {input.fasta} -i chromsizes) -empty | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4+$5}}' > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_{factor}-chipseq-counts-{strand}.bedgraph",
        bam_experimental = "alignment/{sample}_{factor}-chipseq-uniquemappers-experimental.bam",
        bam_spikein = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-spikein.bam" if wc.norm=="spikenorm" and wc.sample in CHIPS else [],
        input_bam_experimental = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-experimental.bam".format(sample=CHIPS[wc.sample]["input"], factor=wc.factor) if wc.norm=="spikenorm" and wc.sample in CHIPS else [],
        input_bam_spikein = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-spikein.bam".format(sample=CHIPS[wc.sample]["input"], factor=wc.factor) if wc.norm=="spikenorm" and wc.sample in CHIPS else []
    output:
        normalized = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus|protection|midpoints"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{strand}-{factor}.log"
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

rule subtract_inputs:
    input:
        ip_sample = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bedgraph",
        input_sample = lambda wc: f"coverage/{wc.norm}/{{sample}}_{FACTOR}-chipseq-{wc.norm}-{wc.strand}.bedgraph".format(sample=CHIPS[wc.sample]["input"]),
    output:
        "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}-input-subtracted.bedgraph",
    log:
        "logs/subtract_inputs/subtract_inputs-{sample}-{norm}-{strand}-{factor}.log"
    shell: """
        (bedtools unionbedg -i {input.ip_sample} {input.input_sample} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4-$5}}' > {output}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bg = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bedgraph",
        fasta = lambda wc: os.path.abspath(config["spike_in"]["fasta"]) if wc.norm=="sicounts" else os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bw"
    params:
        stranded = lambda wc: [] if wc.strand not in ["SENSE", "ANTISENSE"] else """| awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2; print $1"-minus", $2}}' | LC_COLLATE=C sort -k1,1"""
    log :
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{norm}-{strand}-{factor}.log"
    shell: """
        (bedGraphToBigWig {input.bg} <(faidx {input.fasta} -i chromsizes {params.stranded}) {output}) &> {log}
        """

