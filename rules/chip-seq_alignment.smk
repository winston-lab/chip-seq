#!/usr/bin/env python

localrules:
    build_combined_genome,
    bowtie2_build,
    index_bam,

basename = "{exp_name}_{exp_fasta}_{si_name}_{si_fasta}".format(exp_name = config["genome"]["name"],
                                                                exp_fasta = os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0],
                                                                si_name = config["spike_in"]["name"],
                                                                si_fasta = os.path.splitext(os.path.basename(config["spike_in"]["fasta"]))[0]) if SISAMPLES else os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0]

rule build_combined_genome:
    input:
        experimental = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        spikein = config["spike_in"]["fasta"] if SISAMPLES else []
    output:
        "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(build_annotations(config["genome"]["fasta"])))[0], bn=basename),
    params:
        exp_name = config["genome"]["name"],
        si_name = config["spike_in"]["name"] if SISAMPLES else []
    log: "logs/build_combined_genome.log"
    shell: """
        (sed 's/>/>{params.exp_name}_/g' {input.experimental} | \
        cat - <(sed 's/>/>{params.si_name}_/g' {input.spikein}) > {output}) &> {log}
        """

rule bowtie2_build:
    input:
        "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(build_annotations(config["genome"]["fasta"])))[0], bn=basename) if SISAMPLES else os.path.abspath(build_annotations(config["genome"]["fasta"])),
    output:
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["bowtie2"]["index-path"]
    conda: "../envs/bowtie2.yaml"
    log: "logs/bowtie2_build-{basename}.log"
    shell: """
        (bowtie2-build {input} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand("{directory}/{bn}.{num}.bt2", directory = config["bowtie2"]["index-path"],
                                             bn = basename,
                                             num = [1,2,3,4]),
        expand("{directory}/{bn}.rev.{num}.bt2", directory = config["bowtie2"]["index-path"],
                                             bn = basename,
                                             num = [1,2]),
        fastq = f"fastq/cleaned/{{sample}}_{FACTOR}-chipseq-cleaned.fastq.gz",
    output:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
        unaligned_fastq = f"fastq/{{sample}}_{FACTOR}-chipseq-unaligned.fastq.gz",
        log = "logs/align/align_{sample}.log"
    params:
        idx_path = config["bowtie2"]["index-path"],
        minmapq = config["bowtie2"]["minmapq"]
    conda: "../envs/bowtie2.yaml"
    threads : config["threads"]
    shell: """
        (bowtie2 -x {params.idx_path}/{basename} -U {input.fastq} --un-gz {output.unaligned_fastq} -p {threads} | \
         samtools view -buh -q {params.minmapq} - | \
         samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        """

#indexing is required for separating species by samtools view
rule index_bam:
    input:
        f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam"
    output:
        f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input}) &> {log}
        """

rule bam_separate_species:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
        bai = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam.bai",
        fasta = "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(build_annotations(config["genome"]["fasta"])))[0], bn=basename) if SISAMPLES else [],
    output:
        f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-{{species}}.bam",
    params:
        filterprefix = lambda wc: config["spike_in"]["name"] if wc.species=="experimental" else config["genome"]["name"],
        prefix = lambda wc: config["genome"]["name"] if wc.species=="experimental" else config["spike_in"]["name"]
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(faidx {input.fasta} -i chromsizes | \
         grep {params.prefix}_ | \
         awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | \
         grep -v -e 'SN:{params.filterprefix}_' | \
         sed 's/{params.prefix}_//g' | \
         samtools view -bh -@ {threads} -o {output} -) &> {log}
        """

