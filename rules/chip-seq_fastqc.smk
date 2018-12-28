#!/usr/bin/env python

localrules:
    make_barcode_file,
    fastqc_aggregate

# make a file containing adapters for fastQC to search
# barcodes include the 'A' tail
rule make_barcode_file:
    output:
        "fastq/barcodes.tsv"
    run:
        with open(output[0], "w") as out:
            for sample, metadata in SAMPLES.items():
                out.write(f'{sample}\t{metadata["barcode"]}T\n')

rule fastqc_pre_and_unaligned:
    input:
        fastq = lambda wc: {"raw": config["unmatched"] if wc.sample=="unmatched" else SAMPLES[wc.sample]["fastq"],
                            "cleaned": f"fastq/cleaned/{wc.sample}_{FACTOR}-chipseq-cleaned.fastq.gz",
                            "unaligned": f"fastq/{wc.sample}_{FACTOR}-chipseq-unaligned.fastq.gz"
                            }.get(wc.read_status, "KEYERROR"),
        adapters = lambda wc: [] if wc.sample != "unmatched" else "fastq/barcodes.tsv"
    output:
        "qual_ctrl/fastqc/{read_status}/{sample}_fastqc-data-{read_status}.txt"
    params:
        fname = lambda wc: {"raw": "allunmatched" if wc.sample=="unmatched" else re.split('.fq|.fastq', os.path.split(SAMPLES[wc.sample]["fastq"])[1])[0],
                            "cleaned": f"{wc.sample}_{FACTOR}-chipseq-cleaned",
                            "unaligned": f"{wc.sample}_{FACTOR}-chipseq-unaligned"
                    }.get(wc.read_status, "KEYERROR"),
    wildcard_constraints:
        read_status="raw|cleaned|unaligned"
    threads:
        config["threads"]
    log:
        "logs/fastqc/fastqc_{read_status}_{sample}.log"
    run:
        if wildcards.sample=="unmatched":
            shell("""(mkdir -p qual_ctrl/fastqc/raw) &> {log};
                    (cat {input.fastq} > qual_ctrl/fastqc/raw/allunmatched.fastq.gz) &>> {log};
                    (fastqc --adapters {input.adapters} --nogroup --noextract -t {threads} -o qual_ctrl/fastqc/raw qual_ctrl/fastqc/raw/allunmatched.fastq.gz) &>> {log};
                    (unzip -p qual_ctrl/fastqc/raw/{params.fname}_fastqc.zip {params.fname}_fastqc/fastqc_data.txt > {output}) &>> {log}
                    (rm qual_ctrl/fastqc/raw/allunmatched.fastq.gz) &>> {log}""")
        else:
            adapter = SAMPLES[wildcards.sample]["barcode"]+"T"
            shell("""(mkdir -p qual_ctrl/fastqc/{wildcards.read_status}) &> {log};
                    (fastqc --adapters <(echo -e "adapter\t{adapter}") --nogroup --noextract -t {threads} -o qual_ctrl/fastqc/{wildcards.read_status} {input.fastq}) &>> {log};
                    (unzip -p qual_ctrl/fastqc/{wildcards.read_status}/{params.fname}_fastqc.zip {params.fname}_fastqc/fastqc_data.txt > {output}) &>> {log}""")

rule fastqc_aligned:
    input:
        f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam"
    params:
        fname = f"{{sample}}_{FACTOR}-chipseq-uniquemappers",
        adapter = lambda wc: SAMPLES[wc.sample]["barcode"]+"T"
    output:
        "qual_ctrl/fastqc/unique_mappers/{sample}_fastqc-data-unique_mappers.txt"
    threads:
        config["threads"]
    log:
        "logs/fastqc/fastqc_unique_mappers_{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/unique_mappers) &> {log}
        (bedtools bamtofastq -fq qual_ctrl/fastqc/unique_mappers/{params.fname}.fastq -i {input}) &>> {log}
        (fastqc --adapters <(echo -e "adapter\t{params.adapter}") --nogroup --noextract -t {threads} -o qual_ctrl/fastqc/unique_mappers qual_ctrl/fastqc/unique_mappers/{params.fname}.fastq) &>> {log}
        (unzip -p qual_ctrl/fastqc/unique_mappers/{params.fname}_fastqc.zip {params.fname}_fastqc/fastqc_data.txt > {output}) &>> {log}
        (rm qual_ctrl/fastqc/unique_mappers/{params.fname}.fastq) &>> {log}
        """

fastqc_dict = {
        "per_base_qual":
        {   "title" : "Per base sequence quality",
            "fields": "base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus"
            } ,
        "per_tile_qual":
        {   "title" : "Per tile sequence quality",
            "fields": "tile\tbase\tmean\tsample\tstatus"
            },
        "per_seq_qual":
        {   "title" : "Per sequence quality scores",
            "fields": "quality\tcount\tsample\tstatus"
            },
        "per_base_seq_content":
        {   "title" : "Per base sequence content",
            "fields": "base\tg\ta\tt\tc\tsample\tstatus"
            },
        "per_seq_gc":
        {   "title" : "Per sequence GC content",
            "fields": "gc_content\tcount\tsample\tstatus"
            },
        "per_base_n":
        {   "title" : "Per base N content",
            "fields": "base\tn_count\tsample\tstatus"
            },
        "seq_length_dist":
        {   "title" : "Sequence Length Distribution",
            "fields": "length\tcount\tsample\tstatus"
            },
        "seq_duplication":
        {   "title" : "Total Deduplicated Percentage",
            "fields": "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus"
            },
        "adapter_content":
        {   "title" : "Adapter Content",
            "fields": "position\tpct\tsample\tstatus"
            }
        }

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}_fastqc-data-raw.txt", sample=(["unmatched"] if config["unmatched"] else []) + list(SAMPLES.keys())),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}_fastqc-data-cleaned.txt", sample=SAMPLES),
        aligned = expand("qual_ctrl/fastqc/unique_mappers/{sample}_fastqc-data-unique_mappers.txt", sample=SAMPLES),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}_fastqc-data-unaligned.txt", sample=SAMPLES),
    output:
        per_base_qual = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_quality.tsv',
        per_tile_qual = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_tile_quality.tsv',
        per_seq_qual =  f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_sequence_quality.tsv',
        per_base_seq_content = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_sequence_content.tsv',
        per_seq_gc = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_sequence_gc.tsv',
        per_base_n = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_n.tsv',
        seq_length_dist = f'qual_ctrl/fastqc/{FACTOR}-chipseq-sequence_length_distribution.tsv',
        seq_duplication = f'qual_ctrl/fastqc/{FACTOR}-chipseq-sequence_duplication_levels.tsv',
        adapter_content = f'qual_ctrl/fastqc/{FACTOR}-chipseq-adapter_content.tsv',
    run:
        shell("""rm -f {output}""")
        for fastqc_metric, out_path in output.items():
            title = fastqc_dict[fastqc_metric]["title"]
            fields = fastqc_dict[fastqc_metric]["fields"]
            for read_status, read_status_data in input.items():
                for sample_id, fastqc_data in zip((["unmatched"] if config["unmatched"] and read_status=="raw" else []) + list(SAMPLES.keys()), read_status_data):
                    if sample_id=="unmatched" and title=="Adapter Content":
                        shell("""awk 'BEGIN{{FS=OFS="\t"}} /{title}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{m=$2;for(i=2;i<=NF-2;i++)if($i>m)m=$i; print $1, m, "{sample_id}", "{read_status}"}}' {fastqc_data} | tail -n +2 >> {out_path}""")
                    else:
                        shell("""awk 'BEGIN{{FS=OFS="\t"}} /{title}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, "{sample_id}", "{read_status}"}}' {fastqc_data} | tail -n +2 >> {out_path}""")
            shell("""sed -i "1i {fields}" {out_path}""")

rule plot_fastqc_summary:
    input:
        seq_len_dist = f'qual_ctrl/fastqc/{FACTOR}-chipseq-sequence_length_distribution.tsv',
        per_tile = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_tile_quality.tsv',
        per_base_qual = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_quality.tsv',
        per_base_seq = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_sequence_content.tsv',
        per_base_n = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_n.tsv',
        per_seq_gc = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_sequence_gc.tsv',
        per_seq_qual = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_sequence_quality.tsv',
        adapter_content = f'qual_ctrl/fastqc/{FACTOR}-chipseq-adapter_content.tsv',
        seq_dup = f'qual_ctrl/fastqc/{FACTOR}-chipseq-sequence_duplication_levels.tsv',
    output:
        seq_len_dist = f'qual_ctrl/fastqc/{FACTOR}-chipseq-sequence_length_distribution.svg',
        per_tile = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_tile_quality.svg',
        per_base_qual = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_quality.svg',
        per_base_seq = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_sequence_content.svg',
        per_seq_gc = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_sequence_gc.svg',
        per_seq_qual = f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_sequence_quality.svg',
        adapter_content = f'qual_ctrl/fastqc/{FACTOR}-chipseq-adapter_content.svg',
        seq_dup = f'qual_ctrl/fastqc/{FACTOR}-chipseq-sequence_duplication_levels.svg',
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/fastqc_summary.R"

