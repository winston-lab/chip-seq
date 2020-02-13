#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import pyBigWig as pybw

#given paths to bigwig files representing replicates,
#return a dictionary where keys are chromosome names and
#values are the average coverage across replicates
def average_bigwigs(coverage_paths):
    coverage = {}
    for index, path in enumerate(coverage_paths):
        bw = pybw.open(path)
        chroms = bw.chroms()
        for chrom in chroms:
            if index==0:
                coverage[chrom] = bw.values(chrom, 0, chroms[chrom], numpy=True)
            else:
                coverage[chrom] = np.add(bw.values(chrom, 0, chroms[chrom], numpy=True), coverage[chrom])
            if index==len(coverage_paths):
                coverage[chrom] = np.divide(coverage[chrom], index)
        bw.close()
    return coverage

#given "unstranded" chromosome coordinates,
#return 0-based offset of summit position from start.
#If multiple positions have the same max signal, return the mean position
def get_summit(row, coverage):
    local_coverage = coverage[row['chrom']][row['start']:row['end']]
    if not np.any(np.isfinite(local_coverage)):
        return int(len(local_coverage) / 2)
    return int(np.mean(np.argwhere(local_coverage==np.amax(local_coverage[np.isfinite(local_coverage)]))))

def main(condition_paths,
        control_paths,
        diffexp_path,
        narrowpeak_out,
        bed_out):

    #condition and control coverage are imported separately and
    #averaged across replicates in case the number of samples
    #in each group is different
    condition_coverage = average_bigwigs(condition_paths)
    coverage = average_bigwigs(control_paths)
    for chrom in coverage:
        coverage[chrom] = np.add(coverage[chrom], condition_coverage[chrom])

    #we only need to perform operations using start and end as integers,
    #so everything else can be treated as an object to avoid reformatting
    diffexp_df = pd.read_csv(diffexp_path, sep="\t",
                             dtype={'chrom':str,
                                    'start':np.uint32,
                                    'end':np.uint32,
                                    'name':str,
                                    'score':str,
                                    'strand':str,
                                    'log2FC_enrichment':str,
                                    'lfc_SE':str,
                                    'stat':str,
                                    'log10_pval':str,
                                    'log10_padj':str,
                                    'mean_counts':str,
                                    'condition_enrichment':str,
                                    'condition_enrichment_SE':str,
                                    'control_enrichment':str,
                                    'control_enrichment_SE':str})

    if diffexp_df.shape[0] > 0:
        diffexp_df['summit'] = diffexp_df.apply(get_summit, coverage=coverage, axis=1)
        diffexp_df = diffexp_df.assign(summit_start = diffexp_df['start'] + diffexp_df['summit'])
        diffexp_df = diffexp_df.assign(summit_end = diffexp_df['summit_start'] + 1)

    #NOTE: we convert NAs (found in pvalue and score columns) to zero for narrowpeak compatibility
    diffexp_df.to_csv(narrowpeak_out,
                      sep="\t",
                      columns=(['chrom', 'start', 'end', 'name', 'score', 'strand',
                               'log2FC_enrichment', 'log10_pval', 'log10_padj', 'summit'] if
                               diffexp_df.shape[0] > 0 else []),
                      header=False,
                      index=False,
                      float_format="%.3f",
                      encoding='utf-8',
                      na_rep="0")
    diffexp_df.to_csv(bed_out,
                      sep="\t",
                      columns=(['chrom', 'summit_start', 'summit_end', 'name', 'score', 'strand'] if
                          diffexp_df.shape[0] > 0 else []),
                      header=False,
                      index=False,
                      float_format="%.3f",
                      encoding='utf-8',
                      na_rep="0")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Add back summit information to ChIP-seq differential binding results.')
    parser.add_argument('-i', dest = 'condition_paths', type=str, nargs='+', help='BigWigs for all condition samples')
    parser.add_argument('-j', dest = 'control_paths', type=str, nargs='+', help='BigWigs for all control samples')
    parser.add_argument('-d', dest = 'diffexp_path', type=str, help='differential binding results file')
    parser.add_argument('-n', dest = 'narrowpeak_out', type=str, help='output path for narrowPeak file')
    parser.add_argument('-b', dest = 'bed_out', type=str, help='output path for BED file of summit positions')
    args = parser.parse_args()

    main(args.condition_paths,
         args.control_paths,
         args.diffexp_path,
         args.narrowpeak_out,
         args.bed_out)

