#!/usr/bin/env python

localrules:
    combine_peaks

rule callpeaks_macs2:
    input:
        chip_bam = "alignment/{sample}_{factor}-chipseq-uniquemappers-{species}.bam",
        input_bam = lambda wc: expand("alignment/{sample}_{{factor}}-chipseq-uniquemappers-{{species}}.bam", sample=SAMPLES[wc.sample]["control"]),
        fasta = lambda wc: os.path.abspath(build_annotations(config["genome"]["fasta"])) if wc.species=="experimental" else os.path.abspath(config["spike_in"]["fasta"]),
    output:
        tsv = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_peaks.xls",
        peaks = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_peaks.narrowPeak",
        summits = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_summits.bed",
        script = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_model.r",
        pdf = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_model.pdf",
        treat_bg = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_treat_pileup.bdg",
        cntrl_bg = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipseq_control_lambda.bdg"
    params:
        bw = config["peakcalling"]["bw"],
        slocal = config["peakcalling"]["slocal"],
        llocal = config["peakcalling"]["llocal"],
        mfold_low = config["peakcalling"]["mfold_low"],
        mfold_high = config["peakcalling"]["mfold_high"],
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/callpeaks_macs2/callpeaks_macs2_{sample}-{species}-{factor}.log"
    shell: """
        (macs2 callpeak --treatment {input.chip_bam} --control {input.input_bam} --name peakcalling/sample_peaks/{wildcards.sample}_{wildcards.species}-{wildcards.factor}-chipseq --SPMR --bw {params.bw} --mfold {params.mfold_low} {params.mfold_high} --format BAM --gsize $(faidx {input.fasta} -i chromsizes | awk '{{sum += $2}} END {{print sum}}') --qvalue 1 --slocal {params.slocal} --llocal {params.llocal} --keep-dup auto --bdg --call-summits) &> {log}
        (Rscript {output.script}) &>> {log}
        (sed -i -e 's/peakcalling\/sample_peaks\///g' {output.peaks}) &>> {log}
        (sed -i -e 's/peakcalling\/sample_peaks\///g' {output.summits}) &>> {log}
        """

rule idr:
    input:
        #NOTE: for now we take the first two samples since the IDR script only takes two
        #change this if we find a better way to aggregate results
        lambda wc: ["peakcalling/sample_peaks/" + x + "_{species}-{factor}-chipseq_peaks.narrowPeak".format(**wc) for x in get_samples(passing=True, paired=True, groups=wc.group)][0:2]
    output:
        allpeaks = "peakcalling/{group}/{group}_{species}-{factor}-chipseq-idrpeaks-all.tsv",
        filtered = "peakcalling/{group}/{group}_{species}-{factor}-chipseq-idrpeaks-filtered.tsv",
        narrowpeak = "peakcalling/{group}/{group}_{species}-{factor}-chipseq-idrpeaks.narrowPeak",
        summits = "peakcalling/{group}/{group}_{species}-{factor}-chipseq-idrpeaks-summits.bed",
    params:
        idr = int(-125*log2(config["peakcalling"]["idr"]))
    conda:
        "../envs/idr.yaml"
    log:
        "logs/idr/idr-{group}-{species}-{factor}.log"
    shell: """
        (idr -s {input} --input-file-type narrowPeak --rank q.value -o {output.allpeaks} -l {log} --plot --peak-merge-method max) &> {log}
        (awk '$5>{params.idr} || $9=="inf"' {output.allpeaks} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         tee {output.filtered} | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $10}}' | \
         tee {output.narrowpeak} | \
         awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.summits}) &>> {log}
        """

rule combine_peaks:
    input:
        cond = "peakcalling/{condition}/{condition}_{species}-{factor}-chipseq-idrpeaks.narrowPeak",
        ctrl = "peakcalling/{control}/{control}_{species}-{factor}-chipseq-idrpeaks.narrowPeak",
    output:
        "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_{species}-{factor}-peaks.bed"
    log:
        "logs/combine_peaks/combine_peaks-{condition}-v-{control}-{species}-{factor}.log"
    shell: """
        (bedtools multiinter -i {input} | \
         bedtools merge -i stdin | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, ".", 0, "."}}' | \
         sort -k1,1 -k2,2n > {output}) &> {log}
        """

