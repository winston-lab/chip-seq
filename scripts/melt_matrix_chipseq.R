#!/usr/bin/env/ Rscript
library(argparse)
library(tidyverse)
library(magrittr)

parser = ArgumentParser()
parser$add_argument('-i', '--input', type='character')
parser$add_argument('-r', '--refpt', type='character', nargs="+")
parser$add_argument('-g', '--group', type='character', nargs="+")
parser$add_argument('-s', '--sample', type='character', nargs="+")
parser$add_argument('-t', '--sampletype', type='character', nargs=1)
parser$add_argument('-a', '--annotation', type='character', nargs="+")
parser$add_argument('-b', '--binsize', type='integer')
parser$add_argument('-u', '--upstream', type='integer')
parser$add_argument('-o', '--output', type='character')

args = parser$parse_args()

melt = function(inmatrix, refpt, group, sample, sampletype,
                annotation, binsize, upstream, outpath){
    raw = read_tsv(inmatrix, skip=3, col_names=FALSE)
    names(raw) = seq(ncol(raw))

    df = raw %>%
          rownames_to_column(var="index") %>%
          gather(key = variable, value=value, -index, convert=TRUE) %>%
          filter(!is.na(value)) %>%
          transmute(group = group, sample = sample,
                    sampletype = sampletype, annotation = annotation,
                    index = as.integer(index),
                    position = variable,
                    cpm = as.numeric(value))
    if(binsize>1){
        df %<>% mutate(position = (as.numeric(position)*binsize-(upstream+1.5*binsize))/1000)
    } else if (refpt=="TES"){
        df %<>% mutate(position = (as.numeric(position)-(1+upstream))/1000)
    } else {
        df %<>% mutate(position = (as.numeric(position)-(2+upstream))/1000)
    }
    write_tsv(df, path=outpath, col_names=FALSE)
    return(df)
}

melt(inmatrix = args$input,
     refpt = args$refpt,
     group = paste(args$group, collapse=" "),
     sample = paste(args$sample, collapse=" "),
     sampletype = args$sampletype,
     annotation = paste(args$annotation, collapse=" "),
     binsize = args$binsize,
     upstream = args$upstream,
     outpath = args$output)

