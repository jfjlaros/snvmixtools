#!/usr/bin/env python

import argparse
import matplotlib

matplotlib.use("SVG")
from matplotlib import pyplot

import wiggelen

from . import docSplit, version, usage
from . import snvmix_parse

def freqs(snvmix_handle, picture_handle, log_handle):
    """
    Plot the distribution of minor allele frequencies.

    :arg snvmix_handle:
    :type snvmix_handle: stream
    :arg picture_handle:
    :type picture_handle: stream
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
    pyplot.yscale("log")
    pyplot.axis([0, 0.5, 0.1, 100000])
    pyplot.hist(data, bins=100)
    pyplot.savefig("{}".format(picture_handle.name))
#freqs

def snvmix2wig(snvmix_handle, wiggle_handle):
    """
    Convert an SNVMix file to wiggle.
    """
    wiggelen.write(map(lambda x: (x.chromosome, x.position,
        min(x.reference_count, x.alternative_count)),
        snvmix_parse.walker(snvmix_handle)), track=wiggle_handle)
#snvmix2wig

def main():
    """
    Main entry point.
    """
    input_parser = argparse.ArgumentParser(add_help=False)
    input_parser.add_argument("INPUT", type=argparse.FileType('r'),
        help="input file")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])
    parser.add_argument('-v', action="version", version=version(parser.prog))
    subparsers = parser.add_subparsers(dest="subcommand")

    freqs_parser = subparsers.add_parser("freqs", parents=[input_parser],
        description=docSplit(freqs))
    freqs_parser.add_argument("OUTPUT", type=argparse.FileType('w'), nargs=2,
        help="output files")

    snvmix2wig_parser = subparsers.add_parser("snvmix2wig",
        parents=[input_parser], description=docSplit(snvmix2wig))
    snvmix2wig_parser.add_argument("OUTPUT", type=argparse.FileType('w'),
        help="output file")

    args = parser.parse_args()
    
    if args.subcommand == "freqs":
        freqs(args.INPUT, *args.OUTPUT)

    if args.subcommand == "snvmix2wig":
        snvmix2wig(args.INPUT, args.OUTPUT)
#main

if __name__ == "__main__":
    main()
