#!/usr/bin/env python

import argparse
import numpy as np
import pandas as pd
import pyBigWig as pybw
from scipy import stats
# from statsmodels.stats import multitest


def test_shift(row,
               coverage_objects,
               groups,
               control_id):
    positions_control = []
    positions_condition = []
    for coverage_object, group in zip(coverage_objects, groups):
        local_coverage = coverage_object.values(row["chrom"],
                                                row["start"],
                                                row["end"],
                                                numpy=True)
        if row["strand"] == "-":
            local_coverage = np.flip(local_coverage)
        local_coverage = (local_coverage - np.min(local_coverage)) / np.ptp(local_coverage)
        local_coverage = np.cumsum(local_coverage) / np.sum(local_coverage)
        position = np.argmax(local_coverage > 0.5)
        if group == control_id:
            positions_control.append(position)
        else:
            positions_condition.append(position)

    position_mean_control = np.mean(positions_control)
    position_mean_condition = np.mean(positions_condition)
    position_sd_control = np.std(positions_control)
    position_sd_condition = np.std(positions_condition)
    shift_estimate = position_mean_condition - position_mean_control

    tstat, pval = stats.ttest_ind(positions_control,
                                  positions_condition,
                                  equal_var=True)

    return pd.Series([position_mean_control,
                      position_sd_control,
                      position_mean_condition,
                      position_sd_condition,
                      shift_estimate,
                      pval])


def main(coverage_paths,
         groups,
         control_id,
         annotation_path,
         output_path):
    annotations = pd.read_csv(annotation_path,
                              sep="\t",
                              names=["chrom",
                                     "start",
                                     "end",
                                     "name",
                                     "score",
                                     "strand"],
                              usecols=list(range(6)))

    coverage_objects = [pybw.open(x) for x in coverage_paths]

    annotations[["position_mean_control",
                 "position_sd_control",
                 "position_mean_condition",
                 "position_sd_condition",
                 "shift_estimate",
                 "shift_pval"]] = annotations.apply(test_shift,
                                                    coverage_objects=coverage_objects,
                                                    groups=groups,
                                                    control_id=control_id,
                                                    axis=1)
    annotations.to_csv(output_path,
                       sep="\t",
                       index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Get shift results.")
    parser.add_argument('-i', dest='coverage_paths', type=str, nargs='+',
                        help='Path to coverage files.')
    parser.add_argument('-g', dest='groups', type=str, nargs='+',
                        help='Experimental groups of coverage files.')
    parser.add_argument('-c', dest='control_id', type=str,
                        help='Control experimental group.')
    parser.add_argument('-b', dest='annotation_path', type=str,
                        help='Path to BED6+ file of annotations to find shifts over.')
    parser.add_argument('-o', dest='output_path', type=str,
                        help='Path to write shift results to.')
    args = parser.parse_args()
    main(coverage_paths=args.coverage_paths,
         groups=args.groups,
         control_id=args.control_id,
         annotation_path=args.annotation_path,
         output_path=args.output_path)
