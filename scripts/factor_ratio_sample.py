#!/usr/bin/env python

import argparse
import numpy as np
import pyBigWig as pybw
from itertools import compress

def get_weights_from_groups(groups):
    return [1 / groups.count(x) for x in groups]

def get_weighted_global_mean(coverage_paths,
        weights):
    numerator = 0
    denominator = 0
    for index, coverage_path in enumerate(coverage_paths):
        cumulative_signal = 0
        n_finite_bases = 0
        bw = pybw.open(coverage_path)
        for chrom, size in bw.chroms().items():
            coverage = bw.values(chrom, 0, size, numpy=True)
            finite_indices = np.isfinite(coverage)
            n_finite_bases += np.sum(finite_indices)
            cumulative_signal += np.sum(coverage[finite_indices])
        bw.close()
        numerator += weights[index] * cumulative_signal
        denominator += weights[index] * n_finite_bases
    return numerator / denominator

def get_weighted_global_sdev(coverage_paths,
        weights,
        global_mean):
    numerator = 0
    denominator = 0
    for index, coverage_path in enumerate(coverage_paths):
        cumulative_signal = 0
        n_finite_bases = 0
        bw = pybw.open(coverage_path)
        for chrom, size in bw.chroms().items():
            coverage = bw.values(chrom, 0, size, numpy=True)
            finite_indices = np.isfinite(coverage)
            n_finite_bases += np.sum(finite_indices)
            cumulative_signal += np.sum(np.power(np.subtract(coverage[finite_indices], global_mean), 2))
        bw.close()
        numerator += weights[index] * cumulative_signal
        denominator += weights[index] * n_finite_bases
    return np.sqrt(numerator / denominator)

def average_coverage(coverage_paths):
    coverage = {}
    for index, coverage_path in enumerate(coverage_paths):
        bw = pybw.open(coverage_path)
        for chrom, size in bw.chroms().items():
            if index==0:
                coverage[chrom] = bw.values(chrom, 0, size, numpy=True)
            else:
                coverage[chrom] = np.add(bw.values(chrom, 0, size, numpy=True), coverage[chrom])
            if index==len(coverage_paths):
                coverage[chrom] = np.divide(coverage[chrom], index)
        bw.close()
    return coverage

def standardize_coverage(coverage_path,
        global_mean,
        global_sdev):
    standardized_coverage = {}

    bw = pybw.open(coverage_path)
    for chrom, size in bw.chroms().items():
        values = bw.values(chrom, 0, size, numpy=True)
        standardized_coverage[chrom] = np.divide(np.subtract(values, global_mean), global_sdev)

    bw.close()
    return standardized_coverage

def get_standard_differences(
        numerator_dict,
        denominator_dict):
    differences = {}
    for chrom, numerator_values in numerator_dict.items():
        denominator_values = denominator_dict[chrom]
        differences[chrom] = np.subtract(numerator_values, denominator_values)
    return differences

def main(
        numerator_ratio_paths = ["non-depleted-Spn1-IP-1_Spn1-chipseq-spikenorm-ratio.bw",
            "non-depleted-Spn1-IP-2_Spn1-chipseq-spikenorm-ratio.bw",
            "depleted-Spn1-IP-1_Spn1-chipseq-spikenorm-ratio.bw",
            "depleted-Spn1-IP-2_Spn1-chipseq-spikenorm-ratio.bw"],
        numerator_groups = ["non-depleted",
                "non-depleted",
                "depleted",
                "depleted"],
        denominator_ratio_paths = ["non-depleted-Rpb1-IP-1_Rpb1-chipseq-spikenorm-ratio.bw",
            "non-depleted-Rpb1-IP-2_Rpb1-chipseq-spikenorm-ratio.bw",
            "depleted-Rpb1-IP-1_Rpb1-chipseq-spikenorm-ratio.bw",
            "depleted-Rpb1-IP-2_Rpb1-chipseq-spikenorm-ratio.bw"],
        denominator_groups = ["non-depleted",
            "non-depleted",
            "depleted",
            "depleted"],
        numerator_sample_path = "non-depleted-Spn1-IP-1_Spn1-chipseq-spikenorm-ratio.bw",
        denominator_sample_path = "non-depleted-Rpb1-IP-1_Rpb1-chipseq-spikenorm-ratio.bw",
        bigwig_out = "test2.bw"):

    assert len(numerator_groups) == len(numerator_ratio_paths), "Number of numerator files does not match number of numerator group specifications."
    assert len(denominator_groups) == len(denominator_ratio_paths), "Number of denominator files does not match number of denominator group specifications."

    numerator_weights = get_weights_from_groups(numerator_groups)

    denominator_weights = get_weights_from_groups(denominator_groups)

    numerator_mean = get_weighted_global_mean(
            numerator_ratio_paths,
            numerator_weights)

    denominator_mean = get_weighted_global_mean(
            denominator_ratio_paths,
            denominator_weights)

    numerator_sdev = get_weighted_global_sdev(
            numerator_ratio_paths,
            numerator_weights,
            numerator_mean)

    denominator_sdev = get_weighted_global_sdev(
            denominator_ratio_paths,
            denominator_weights,
            denominator_mean)

    standard_numerator_coverage = standardize_coverage(
            numerator_sample_path,
            numerator_mean,
            numerator_sdev)

    standard_denominator_coverage = standardize_coverage(
            denominator_sample_path,
            denominator_mean,
            denominator_sdev)

    differences = get_standard_differences(
            standard_numerator_coverage,
            standard_denominator_coverage)

    bw_template = pybw.open(numerator_ratio_paths[0])
    bw_out = pybw.open(bigwig_out, "w")
    bw_out.addHeader(list(bw_template.chroms().items()))
    bw_template.close()

    for chrom, values in differences.items():
        bw_out.addEntries(chrom, 0, values=values, span=1, step=1)
    bw_out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get difference in standardized IP/input signal between numerator and denominator factors.')
    parser.add_argument('-n', dest = 'numerator_ratio_paths', type=str, nargs='+', help='Paths to all numerator IP/input ratio coverage bigWigs.')
    parser.add_argument('-d', dest = 'denominator_ratio_paths', type=str, nargs='+', help='Paths to all denominator IP/input ratio coverage bigWigs.')
    parser.add_argument('-j', dest = 'numerator_groups', type=str, nargs='+', help='Experimental groups of numerator coverage files.')
    parser.add_argument('-k', dest = 'denominator_groups', type=str, nargs='+', help='Experimental groups of denominator coverage files.')
    parser.add_argument('-a', "--numerator_sample", dest="numerator_sample_path", type=str, help="Path to numerator IP/input ratio coverage bigWig to normalize.")
    parser.add_argument('-b', "--denominator_sample", dest="denominator_sample_path", type=str, help="Path to denominator IP/input ratio coverage bigWig to normalize.")
    parser.add_argument('-o', dest = 'bigwig_out', type=str, help='Path to output bigWig.')
    args = parser.parse_args()

    main(numerator_ratio_paths = args.numerator_ratio_paths,
         numerator_groups = args.numerator_groups,
         denominator_ratio_paths = args.denominator_ratio_paths,
         denominator_groups = args.denominator_groups,
         numerator_sample_path = args.numerator_sample_path,
         denominator_sample_path = args.denominator_sample_path,
         bigwig_out = args.bigwig_out)

