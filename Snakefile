#!/usr/bin/env python

import os
import re
import itertools
from math import log2

configfile: "config.yaml"

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

configfile: build_annotations("config.yaml")

FACTOR = config["factor"]

INPUTS = config["input_samples"]
INPUTS_PASSING = {k:v for k,v in INPUTS.items() if v["pass-qc"]}
INPUTS_SISAMPLES = {k:v for k,v in INPUTS.items() if v["spikein"]}
INPUTS_SIPASSING = {k:v for k,v in INPUTS_SISAMPLES.items() if v["pass-qc"]}
CHIPS = config["chip_samples"]
CHIPS_PASSING = {k:v for k,v in CHIPS.items() if v["pass-qc"]}
CHIPS_SISAMPLES = {k:v for k,v in CHIPS.items() if v["spikein"]}
CHIPS_SIPASSING = {k:v for k,v in CHIPS_SISAMPLES.items() if v["pass-qc"]}
SAMPLES = {**INPUTS, **CHIPS}
PASSING = {**INPUTS_PASSING, **CHIPS_PASSING}
SISAMPLES = {**INPUTS_SISAMPLES, **CHIPS_SISAMPLES}
SIPASSING = {**INPUTS_SIPASSING, **CHIPS_SIPASSING}
GROUPS = set(v["group"] for (k,v) in CHIPS.items())

#groups which have at least one passing chip and input sample, so that they are valid for peakcalling and differential binding
validgroups = set(v["group"] for k,v in CHIPS_PASSING.items() if v["input"] in INPUTS_PASSING)
validgroups_si = set(v["group"] for k,v in CHIPS_SIPASSING.items() if v["input"] in INPUTS_SIPASSING)

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))

comparisons_si =  config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups_si]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups_si]))

FIGURES = config["figures"]

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys()) + ["unmatched"]),
    group = "|".join(set(re.escape(v["group"]) for k,v in CHIPS.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + (conditiongroups_si if comparisons_si else []) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + (controlgroups_si if comparisons_si else []) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys())),
    annotation = "|".join(re.escape(x) for x in set(list(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES])) + list(config["differential_occupancy"]["annotations"].keys() if config["differential_occupancy"]["annotations"] else []) + ["peaks"])),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    strand = "|".join(list(itertools.chain.from_iterable([[x, x+"-input-subtracted"] for x in ["plus", "minus", "midpoints", "protection"]]))),
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    factor=FACTOR

status_norm_sample_dict = {
    "all":
        {   "libsizenorm" : SAMPLES,
            "spikenorm" : SISAMPLES
        },
    "passing":
        {   "libsizenorm" : PASSING,
            "spikenorm" : SIPASSING
        }
    }

def get_samples(status, norm, groups):
    if "all" in groups:
        return(list(status_norm_sample_dict[status][norm].keys()))
    else:
        return([k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in groups])

include: "rules/chip-seq_clean_reads.smk"
include: "rules/chip-seq_alignment.smk"
include: "rules/chip-seq_fastqc.smk"
include: "rules/chip-seq_library_processing_summary.smk"
include: "rules/chip-seq_peakcalling.smk"
include: "rules/chip-seq_genome_coverage.smk"
include: "rules/chip-seq_sample_similarity.smk"
include: "rules/chip-seq_datavis.smk"
include: "rules/chip-seq_differential_binding.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #fastqc
        f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_sequence_content.svg',
        #alignment
        expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bw", sample=SAMPLES, factor=FACTOR, norm=["counts","libsizenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bw", sample=SISAMPLES, factor=FACTOR, norm=["sicounts", "spikenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/libsizenorm/{sample}_{factor}-chipseq-libsizenorm-{strand}.bw", sample=CHIPS, factor=FACTOR, strand=[x+"-input-subtracted" for x in ["plus","minus","protection","midpoints"]]),
        expand("coverage/spikenorm/{sample}_{factor}-chipseq-spikenorm-{strand}.bw", sample=CHIPS_SISAMPLES, factor=FACTOR, strand=[x+"-input-subtracted" for x in ["plus","minus","protection","midpoints"]]),
        f"qual_ctrl/read_processing/{FACTOR}-chipseq_read-processing-loss.svg",
        expand("qual_ctrl/spikein/{factor}-chipseq_spikein-plots-{status}.svg", factor=FACTOR, status=statuscheck(SISAMPLES, SIPASSING)) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}-chipseq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), status=statuscheck(SISAMPLES, SIPASSING), windowsize=config["scatterplot_binsizes"], factor=FACTOR) if SISAMPLES and comparisons_si else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}-chipseq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), status=statuscheck(SAMPLES, PASSING), windowsize=config["scatterplot_binsizes"], factor=FACTOR),
        #datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{factor}}-chipseq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, status=statuscheck(SAMPLES, PASSING), readtype=["midpoints", "midpoints-input-subtracted", "protection", "protection-input-subtracted"], factor=FACTOR) if config["plot_figures"] else [],
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{factor}}-chipseq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), figure=FIGURES, status=statuscheck(SISAMPLES, SIPASSING), readtype=["midpoints", "midpoints-input-subtracted", "protection", "protection-input-subtracted"], factor=FACTOR) if comparisons_si and config["plot_figures"] else [],
        #differential binding
        expand(expand("diff_binding/{{annotation}}/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{factor}}-chipseq-libsizenorm-{{annotation}}-diffbind-results-all.tsv", zip, condition=conditiongroups, control=controlgroups), annotation=list(config["differential_occupancy"]["annotations"].keys() if config["differential_occupancy"]["annotations"] else [])+["peaks"], factor=FACTOR),
        expand(expand("diff_binding/{{annotation}}/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{factor}}-chipseq-spikenorm-{{annotation}}-diffbind-results-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si), annotation=list(config["differential_occupancy"]["annotations"].keys() if config["differential_occupancy"]["annotations"] else [])+["peaks"], factor=FACTOR) if comparisons_si else []

