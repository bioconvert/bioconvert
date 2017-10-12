# -*- coding: utf-8 -*-
#
#  This file is part of Biokit software
#
#  Copyright (c) 2016 - Biokit Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
""".. rubric:: Standalone application dedicated to coverage"""
import os
import shutil
import glob
import sys
import argparse
from optparse import OptionParser
from argparse import RawTextHelpFormatter
import importlib

from bioconvert import logger, bioconvert_debug_level
from bioconvert.converters.registry import Registry


class Options(argparse.ArgumentParser):
    def  __init__(self, prog="convertor"):
        usage = """\nUSAGE

        convertor inputfile.bam outputfile.bam
        convertor inputfile.sam outputfile.bam
        convertor inputfile.fastq outputfile.fasta

        """

        epilog = """ ----    """

        description = """DESCRIPTION:

Convertor infer the formats from the extension name. We do not scan the
input file. Therefore users must ensure that their input format files are
properly formatted.

        """
        super(Options, self).__init__(usage=usage, prog=prog,
                description=description, epilog=epilog)

        # options to fill the config file
        group = self.add_argument_group("Optional arguments")
        group.add_argument("-f", "--formats", dest="format",
            action="store_true", default=False,
            help=("List available format ."))
        group.add_argument("-l", "--logging-level", dest="logging_level",
            default="INFO",
            help=("List available format ."))
        group.add_argument("-x", "--input-format", dest="input_format",
            default=None,
            help=("Provide the input format. Check the --formats to see valid input name"))


def main(args=None):
    from easydev.console import purple, underline
    print(purple("Welcome to bioconvert (bioconvert.readthedocs.io)"))
    mapper = Registry()

    if args is None:
        args = sys.argv[:]

    user_options = Options(prog="converter")

    # If --help or no options provided, show the help
    if "-f" in args or "--formats" in args:
        options = user_options.parse_args(args[1:])
        if options.format:
            print("Available mapping:")
            print("==================")
            for k in sorted(mapper.get_conversions()):
                print("{} -> {}".format(k[0], k[1]))
            sys.exit(0)

    if len(args) < 3:
        user_options.parse_args(["prog", "--help"])
    else:
        infile = args[1]
        outfile = args[2]
        options = user_options.parse_args(args[3:])

    # Set the logging level
    bioconvert_debug_level(options.logging_level)

    # Users may provide information about the input file.
    # Indeed, the input may be a FastQ file but with an extension
    # that is not standard. For instance fq instead of fastq
    # If so, we can use the --input-format fastq to overwrite the
    # provided filename extension
    inext = "." + os.path.splitext(infile)[-1][1:]
    outext = "." + os.path.splitext(outfile)[-1][1:]

    if options.input_format:
        inext = options.input_format
        if inext[0] != ".":
            inext = "." + inext

    # From the input parameters 1 and 2, we get the module name
    try:
        logger.info("Input: {}".format(inext))
        logger.info("Output: {}".format(outext))
        class_converter = mapper[(inext, outext)]
    except KeyError:
        print(mapper)
        print(inext)
        print(outext)

        # Is the module name available in biokit ? If not, let us tell the user
        msg = "Request input format ({}) to output format (({}) is not available in converters"
        logger.critical(msg.format(inext, outext))
        logger.critical("Use --formats to know the available formats")
        sys.exit(1)


    # If the module exists, it is part of the MapperRegitry dictionary and
    # we should be able to import it dynamically, create the class and call
    # the instance
    logger.info("Converting from {} to {}".format(inext, outext))
    convert = class_converter(infile, outfile)
    convert()
    logger.info("Done")


if __name__ == "__main__":
   main()

