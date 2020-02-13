#!/usr/bin/env python

import argparse
import numpy as np
import pyBigWig as pybw

def main(chip_in="non-depleted-Rpb1-IP-3_Rpb1-chipseq-spikenorm-midpoints_smoothed.bw",
         input_in="non-depleted-untagged-input-3_Rpb1-chipseq-spikenorm-midpoints_smoothed.bw",
         ratio_out="ratio.bw"):
    chip = pybw.open(chip_in)
    input = pybw.open(input_in)
    ratio = pybw.open(ratio_out, "w")

    assert chip.chroms() == input.chroms(), "ChIP and input bigWig chromosomes don't match."

    ratio.addHeader(list(chip.chroms().items()))

    for chrom in chip.chroms():
        chip_values = chip.values(chrom, 0, chip.chroms(chrom), numpy=True)
        input_values = input.values(chrom, 0, chip.chroms(chrom), numpy=True)
        ratio.addEntries(chrom, 0, values=np.log2(np.divide(chip_values, input_values)), span=1, step=1)

    chip.close()
    input.close()
    ratio.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Given two bigWig coverage files, generate a coverage file of their log2 ratio.')
    parser.add_argument('-c', dest='chip_in', type=str, help='Path to numerator (ChIP) bigWig.')
    parser.add_argument('-i', dest='input_in', type=str, help='Path to denominator (input) bigWig.')
    parser.add_argument('-o', dest='ratio_out', type=str, help='Path to output bigWig.')
    args = parser.parse_args()
    main(chip_in=args.chip_in,
         input_in=args.input_in,
         ratio_out=args.ratio_out)

