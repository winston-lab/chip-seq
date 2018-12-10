#!/usr/bin/env python

#   search for and remove barcode sequences at 3' end of read
#   do quality trimming with --nextseq-trim 2-color quality trimming
rule clean_reads:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"],
    output:
        fastq = f"fastq/cleaned/{{sample}}_{FACTOR}-chipseq-cleaned.fastq.gz",
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        adapter = lambda wc: SAMPLES[wc.sample]["barcode"] + "T",
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
    conda:
        "../envs/cutadapt.yaml"
    threads:
        config["threads"]
    shell: """
        (cutadapt --cores={threads} --adapter={params.adapter} --error-rate=0.1 --no-indels --nextseq-trim={params.qual_cutoff} --trim-n --minimum-length=5 --output={output.fastq} {input}) &> {output.log}
        """

