#!/usr/bin/env python

from __future__ import division

import argparse
import matplotlib
import os

from math import *
from matplotlib import pyplot
from scipy import stats

import wiggelen
import pybedtools

from . import ProtectedFileType, doc_split, usage, version
from . import snvmix_parse

plotter = {
    "min": lambda x, y: min(x, y),
    "min_norm": lambda x, y: min(x, y) / (x + y)
}

def freqs(snvmix_handle, output_handle, log_handle, threshold=0,
        filter_function=""):
    """
    Plot the distribution of minor allele frequencies.

    :arg snvmix_handle: Open readable handle to an SNVMix file.
    :type snvmix_handle: stream
    :arg output_handle: Open writable handle to an image file.
    :type output_handle: stream
    :arg log_handle: Open writable handle to a text file.
    :type log_handle: stream
    :arg threshold: Simple coverage threshold.
    :type threshold: int
    :arg filter_function: Custom filter function (overrides *threshold*).
    :type filter_function: str
    """
    filter_func = lambda x, y: x + y > threshold
    if filter_function:
        filter_func = eval("lambda " + filter_function)

    data = []
    total = 0
    for record in snvmix_parse.walker(snvmix_handle):
        if filter_func(record.reference_count, record.alternative_count):
            variants = min(record.reference_count, record.alternative_count)
            data.append(variants / float(record.reference_count +
                record.alternative_count))
            total += variants
        #if
    #for

    log_handle.write("{}\n".format(total))
    pyplot.axis([0, 0.5, 0.1, 100000])
    if data:
        pyplot.hist(data, bins=100, log=True)
    pyplot.savefig("{}".format(output_handle.name))
#freqs

def varcall(snvmix_handle, output_handle, p_value=0.05, error_rate=0.01,
        multiple_testing_correction=False):
    """
    Call variants based on a one-tailed binomial test.

    :arg snvmix_handle: Open readable handle to an SNVMix file.
    :type snvmix_handle: stream
    :arg output_handle: Open writable handle to an SNVMix file.
    :type output_handle: stream
    :arg p_value: p-value threshold.
    :type p_value: float
    :arg error_rate: Base calling error rate.
    :type error_rate: float
    :arg mutliple_testing_correction: Do a Bonferroni correction.
    :type mutliple_testing_correction: bool
    """
    correction = 1
    if multiple_testing_correction:
        correction = sum(1 for _ in snvmix_parse.walker(snvmix_handle))
        snvmix_handle.seek(0)
    #if

    for record in snvmix_parse.walker(snvmix_handle):
        if stats.binom.sf(min(record.reference_count,
                record.alternative_count), record.reference_count +
                record.alternative_count, error_rate) < p_value / correction:
            output_handle.write(str(record))
#varcall

def snvmix2wig(snvmix_handle, output_handle, plot_choice="min",
        plot_function="", threshold=0, filter_function=""):
    """
    Convert an SNVMix file to wiggle.

    :arg snvmix_handle: Open readable handle to an SNVMix file.
    :type snvmix_handle: stream
    :arg output_handle: Open writable handle to an wiggle file.
    :type output_handle: stream
    :arg plot_choice: Select a predefined plotting function.
    :type plot_choice: str
    :arg plot_function: Custom plot function (overrides *plot_choice*).
    :type plot_function: str
    :arg threshold: Simple coverage threshold.
    :type threshold: int
    :arg filter_function: Custom filter function (overrides *threshold*).
    :type filter_function: str
    """
    plot_func = plotter[plot_choice]
    if plot_function:
        plot_func = eval("lambda " + plot_function)

    filter_func = lambda x, y: x + y > threshold
    if filter_function:
        filter_func = eval("lambda " + filter_function)

    wiggelen.write(map(lambda x: (x.chromosome, x.position,
        float(filter_func(x.reference_count, x.alternative_count) and
        plot_func(x.reference_count, x.alternative_count))),
        snvmix_parse.walker(snvmix_handle)), track=output_handle,
        name=os.path.splitext(os.path.basename(output_handle.name))[0])
#snvmix2wig

def intersect(snvmix_handle, bed_handle, output_handle):
    """
    Intersect an SNVMix file with a BED track.

    :arg snvmix_handle: Open readable handle to an SNVMix file.
    :type snvmix_handle: stream
    :arg bed_handle: Open readable handle to a BED file.
    :type bed_handle: stream
    :arg output_handle: Open writable handle to an wiggle file.
    :type output_handle: stream
    """
    walker = snvmix_parse.walker(snvmix_handle)
    bed_track = pybedtools.BedTool(bed_handle.name)

    for bed_record in bed_track:
        snvmix_record = walker.next()
        while (snvmix_record.chromosome != bed_record.chrom or
                snvmix_record.position < bed_record.start):
            snvmix_record = walker.next()
        while (snvmix_record.chromosome == bed_record.chrom and
                snvmix_record.position <= bed_record.end and
                snvmix_record.position >= bed_record.start):
            output_handle.write(str(snvmix_record))
            snvmix_record = walker.next()
        #while
    #for
#intersect

def main():
    """
    Main entry point.
    """
    snvmix_parser = argparse.ArgumentParser(add_help=False)
    snvmix_parser.add_argument("snvmix_handle", metavar="INPUT",
        type=argparse.FileType('r'), help="input snvmix file")

    bed_parser = argparse.ArgumentParser(add_help=False)
    bed_parser.add_argument("bed_handle", metavar="BED",
        type=argparse.FileType('r'), help="BED file")

    output_parser = argparse.ArgumentParser(add_help=False)
    output_parser.add_argument("output_handle", metavar="OUTPUT",
        type=ProtectedFileType('w'), help="output file")

    filter_parser = argparse.ArgumentParser(add_help=False)
    filter_parser.add_argument("-t", dest="threshold", type=int, default=0,
        help='filter treshold (%(type)s default="%(default)s")')
    filter_parser.add_argument("--filter-function", dest="filter_function",
        type=str, default="",
        help='custom filter function (%(type)s default="%(default)s")')

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers(dest="subcommand")

    freqs_parser = subparsers.add_parser("freqs", parents=[snvmix_parser,
        output_parser, filter_parser], description=doc_split(freqs))
    freqs_parser.add_argument("log_handle", metavar="LOG",
        type=ProtectedFileType('w'), help="log file")
    freqs_parser.set_defaults(func=freqs)

    snvmix2wig_parser = subparsers.add_parser("snvmix2wig",
        parents=[snvmix_parser, output_parser, filter_parser],
        description=doc_split(snvmix2wig))
    snvmix2wig_parser.add_argument("-p", dest="plot_choice", type=str,
        choices=plotter, default="min",
        help='plotting function (%(type)s default="%(default)s")')
    snvmix2wig_parser.add_argument("--plot-function", dest="plot_function",
        type=str, default="",
        help='custom plot function (%(type)s default="%(default)s")')
    snvmix2wig_parser.set_defaults(func=snvmix2wig)

    varcall_parser = subparsers.add_parser("varcall", parents=[snvmix_parser,
        output_parser], description=doc_split(varcall))
    varcall_parser.add_argument("-p", dest="p_value", type=float, default=0.05,
        help='p-value threshold (%(type)s default="%(default)s")')
    varcall_parser.add_argument("-e", dest="error_rate", type=float,
        default=0.01,
        help='base calling error rate (%(type)s default="%(default)s")')
    varcall_parser.add_argument("-m", dest="multiple_testing_correction",
        default=False, action="store_true", help='do a Bonferroni correction')
    varcall_parser.set_defaults(func=varcall)

    intersect_parser = subparsers.add_parser("intersect",
        parents=[snvmix_parser, bed_parser, output_parser],
        description=doc_split(intersect))
    intersect_parser.set_defaults(func=intersect)

    try:
        args = parser.parse_args()
    except IOError, error:
        parser.error(error)
    
    try:
        args.func(**{k: v for k, v in vars(args).items()
            if k not in ("func", "subcommand")})
    except ValueError, error:
        parser.error(error)
#main

if __name__ == "__main__":
    main()
