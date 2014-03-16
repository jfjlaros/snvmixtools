#!/usr/bin/env python

import argparse
import matplotlib

matplotlib.use("SVG")
from matplotlib import pyplot

from . import snvmix_parse

def freqs(snvmix_handle, picture_handle, log_handle):
    """
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

def main():
    """
    Main entry point.
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )#description=usage[0], epilog=usage[1])
    #parser.add_argument('-v', action="version", version=version(parser.prog))

    parser.add_argument("INPUT", type=argparse.FileType('r'),
        help="input file")
    parser.add_argument("OUTPUT", type=argparse.FileType('w'), nargs=2,
        help="output files")

    args = parser.parse_args()
    
    freqs(args.INPUT, *args.OUTPUT)
#main

if __name__ == "__main__":
    main()
