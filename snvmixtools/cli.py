#!/usr/bin/env python

from __future__ import division

import argparse
import matplotlib

from matplotlib import pyplot

import wiggelen
import pybedtools

from . import docSplit, version, usage
from . import snvmix_parse

plotter = {
    "min": lambda x, y: min(x, y),
    "min_norm": lambda x, y: min(x, y) / (x + y)
}

def freqs(snvmix_handle, output_handle, log_handle):
    """
    Plot the distribution of minor allele frequencies.

    :arg snvmix_handle:
    :type snvmix_handle: stream
    :arg output_handle:
    :type output_handle: stream
    :arg log_handle:
    :type log_handle: stream
    """
    data = []
    total = 0
    for record in snvmix_parse.walker(snvmix_handle):
        variants = min(record.reference_count, record.alternative_count)
        data.append(variants / float(record.reference_count +
            record.alternative_count))
        total += variants
    #for

    log_handle.write("{}\n".format(total))
    pyplot.axis([0, 0.5, 0.1, 100000])
    pyplot.hist(data, bins=100, log=True)
    pyplot.savefig("{}".format(output_handle.name))
#freqs

def snvmix2wig(snvmix_handle, output_handle, plot_function="min"):
    """
    Convert an SNVMix file to wiggle.
    """
    plot_func = plotter[plot_function]

    wiggelen.write(map(lambda x: (x.chromosome, x.position,
        plot_func(x.reference_count, x.alternative_count)),
        snvmix_parse.walker(snvmix_handle)), track=output_handle)
#snvmix2wig

def intersect(snvmix_handle, bed_handle, output_handle):
    """
    Intersect an SNVMix file with a BED track.
    """
    bed_track = pybedtools.BedTool(bed_handle.name)

    for record in bed_track:
        print record
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
        type=argparse.FileType('w'), help="output file")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers(dest="subcommand")

    freqs_parser = subparsers.add_parser("freqs", parents=[snvmix_parser,
        output_parser], description=docSplit(freqs))
    freqs_parser.add_argument("log_handle", metavar="LOG",
        type=argparse.FileType('w'), help="log file")
    freqs_parser.set_defaults(func=freqs)

    snvmix2wig_parser = subparsers.add_parser("snvmix2wig",
        parents=[snvmix_parser, output_parser],
        description=docSplit(snvmix2wig))
    snvmix2wig_parser.add_argument("-p", dest="plot_function",
        metavar="FUNCTION", choices=plotter, default="min",
        help="plotting function")
    snvmix2wig_parser.set_defaults(func=snvmix2wig)

    intersect_parser = subparsers.add_parser("intersect",
        parents=[snvmix_parser, bed_parser, output_parser],
        description=docSplit(intersect))
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
