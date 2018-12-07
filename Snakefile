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

#groups which have at least one passing chip and input sample, so that they are valid for peakcalling
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
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    strand = "SENSE|ANTISENSE|plus|minus|midpoints|protection",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    factor=FACTOR

# status_norm_sample_dict = {
#     "all":
#         {   "libsizenorm" : SAMPLES,
#             "spikenorm" : SISAMPLES
#         },
#     "passing":
#         {   "libsizenorm" : PASSING,
#             "spikenorm" : SIPASSING
#         }
#     }

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
# include: "rules/mnase-seq_sample_similarity.smk"
# include: "rules/mnase-seq_datavis.smk"
# include: "rules/mnase-seq_quantification.smk"
# include: "rules/mnase-seq_differential_occupancy.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #fastqc
        f'qual_ctrl/fastqc/{FACTOR}-chipseq-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_{factor}-chipseq-uniquemappers.bam", sample=SAMPLES, factor=FACTOR),
        #coverage
        expand("coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bw", sample=SAMPLES, factor=FACTOR, norm=["counts","libsizenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{strand}.bw", sample=SISAMPLES, factor=FACTOR, norm=["sicounts", "spikenorm"], strand=["plus","minus","protection","midpoints"]),
        f"qual_ctrl/read_processing/{FACTOR}-chipseq_read-processing-loss.svg",
        expand("qual_ctrl/spikein/{factor}-chipseq_spikein-plots-{status}.svg", factor=FACTOR, status=["all","passing"]) if SISAMPLES else [],
        #expand("coverage/libsizenorm/{sample}_mnase-midpoint_smoothed-libsizenorm.bw", sample=SAMPLES),
        #expand("coverage/spikenorm/{sample}_mnase-midpoint_smoothed-spikenorm.bw", sample=SISAMPLES),
        ##quality controls
        #"qual_ctrl/read_processing/mnase-seq_read_processing-loss.svg",
        #"qual_ctrl/fragment_length_distributions/mnase-seq_fragment_length_distributions.svg",
        #expand("qual_ctrl/spikein/mnase-seq_spikein-plots-{status}.svg", status=["all","passing"]) if SISAMPLES else [],
        #expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_mnase-seq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all","passing"], windowsize=config["scatterplot_binsizes"]) if SISAMPLES and comparisons_si else [],
        #expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_mnase-seq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], windowsize=config["scatterplot_binsizes"]),
        ##datavis
        #expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-seq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) if config["plot_figures"] and SISAMPLES and comparisons_si else [],
        #expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/mnase-seq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, readtype=["midpoint","wholefrag"], status=["all","passing"]) if config["plot_figures"] else [],
        ##call nucleosomes
        #expand("nucleosome_quantification/{condition}-v-{control}/spikenorm/reference_positions.xls", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING and comparisons_si else [],
        #expand("nucleosome_quantification/{condition}-v-{control}/libsizenorm/reference_positions.xls", zip, condition=conditiongroups, control=controlgroups),
        #expand("nucleosome_quantification/{condition}-v-{control}/spikenorm/{condition}-v-{control}_spikenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING and comparisons_si else [],
        #expand("nucleosome_quantification/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_libsizenorm-dyad-shift-histogram.svg", zip, condition=conditiongroups, control=controlgroups),
        ##danpos over annotations
        #expand(expand("nucleosome_quantification/regions/{{figure}}/spikenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_spikenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups_si, control=controlgroups_si), figure=QUANT) if SIPASSING and comparisons_si else [],
        #expand(expand("nucleosome_quantification/regions/{{figure}}/libsizenorm/{condition}-v-{control}/{{figure}}_{condition}-v-{control}_libsizenorm-individual-occupancy-heatmaps.svg", zip, condition=conditiongroups, control=controlgroups), figure=QUANT),
        ##differential nucleosome levels over transcripts
        #expand("diff_levels/{condition}-v-{control}/spikenorm/{condition}-v-{control}-mnase-seq-results-spikenorm-all.tsv", zip, condition=conditiongroups_si, control=controlgroups_si) if SIPASSING and comparisons_si else [],
        #expand("diff_levels/{condition}-v-{control}/libsizenorm/{condition}-v-{control}-mnase-seq-results-libsizenorm-all.tsv", zip, condition=conditiongroups, control=controlgroups),

