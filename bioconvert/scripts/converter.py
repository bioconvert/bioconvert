# -*- coding: utf-8 -*-

##############################################################################
#  This file is part of Biokit software
#
#  Copyright (c) 2017 - Biokit Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
""".. rubric:: Standalone application dedicated to conversion"""
import os
import sys
import argparse
import colorlog

import bioconvert
from bioconvert.core.registry import Registry
from bioconvert.core.converter import Bioconvert
from bioconvert import extensions



class ConvAction(argparse.Action):

    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help="show all formats available and exit"):
        super(ConvAction, self).__init__(option_strings=option_strings,
                                         dest=dest,
                                         default=default,
                                         nargs=0,
                                         help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        # the -v --verbosity options may not be parsed yet (if located after -f on command line)
        # So I do it myself
        v_nb = ''.join([opt for opt in sys.argv if opt.startswith("-v")]).count('v')

        mapper = Registry()
        print("Available mapping:")
        print("==================")
        for k in sorted(mapper.get_conversions()):
            print("{} -> {}".format(k[0], k[1]))
        sys.exit(0)

def main(args=None):

    if args is None:
        args = sys.argv[1:]

    if "--version" in args:
        print("Bioconvert version {}".format(bioconvert.version))
        sys.exit(0)


    from easydev.console import purple, underline
    if "-v" in args or "--verbosity" in args:
        print(purple("Welcome to bioconvert (bioconvert.readthedocs.io)"))

    arg_parser = argparse.ArgumentParser(prog="bioconvert",
                                         epilog=" ----    ",
                                         description="""Convertor infer the
                                         formats from the extension name. We do
                                         not scan the input file. Therefore
                                         users must ensure that their input
                                         format files are properly
                                         formatted.""",
                                         usage="""
    # convert fastq to fasta
    bioconvert test.fastq test.fasta

    # if input extension is not standard, use -i to specify it
    bioconvert test.FASTQ test.fasta -i fastq

    bioconvert test.fastq -o fasta

    # You may have several inputs, in which case wildcards are possible
    # Note, however, the quotes that are required
    bioconvert "test*.fastq" -o fasta

    # batch is also possible. 
    bioconvert "test*.fastq" -o fasta -m 

    Note the difference between the two previous commands !!


    For more information, please type:

        bioconvert --help
""")
    arg_parser.add_argument("input_file",
            default=None,
            help="The path to the file to convert.")
    arg_parser.add_argument("output_file", nargs="?",
            default=None,
            help="The path where the result will be stored.")

    arg_parser.add_argument("-F", "--formats",
                            action=ConvAction,
                            default=False,
                            help="Display available formats and exit.")
    arg_parser.add_argument("-v", "--verbosity",
                            default="INFO",
                            help="Set the outpout verbosity. Should be one of DEBUG, INFO, WARNING, ERROR, CRITICAL")
    arg_parser.add_argument("--raise-exception",
                            action="store_true",
                            help="Let exception ending the execution be raised and displayed")
    arg_parser.add_argument("-l", "--level", dest="verbosity", 
                            default="INFO",
                            help="same as --verbosity")
    arg_parser.add_argument("-i", "--input-format",
                            default=None,
                            help="Provide the input format. Check the --formats to see valid input name")
    arg_parser.add_argument("-o", "--output-format",
                            default=None,
                            help="Provide the output format. Check the --formats to see valid input name")
    arg_parser.add_argument("-x", "--threads",
                            default=None,
                            type=int,
                            help="Number of threads. Depends on the underlying tool")
    arg_parser.add_argument("-m", "--batch",
                            default=False, action="store_true",
                            help="for batch effect")

    arg_parser.add_argument("-c", "--method",
                            default=None,
                            help="A converter may have several methods")

    arg_parser.add_argument("-f", "--force",
                            action="store_true",
                            help="if outfile exists, it is overwritten with this option")

    arg_parser.add_argument("-s", "--show-methods",
                            default=False,
                            action="store_true",
                            help="A converter may have several methods")

    arg_parser.add_argument("-b", "--benchmark",
                            default=False,
                            action="store_true",
                            help="Running all available methods")

    arg_parser.add_argument("-N", "--benchmark-N",
                            default=5,
                            type=int,
                            help="Number of trials for each methods")

    args = arg_parser.parse_args(args)

    # Set the logging level
    bioconvert.logger.level = args.verbosity


    # Figure out whether we have several input files or not
    # Are we in batch mode ?
    import glob
    if args.batch:
        filenames = glob.glob(args.input_file)
    else:
        filenames = [args.input_file]

    for filename in filenames:
        args.input_file = filename
        try:
            analysis(args)
        except Exception as e:
            if args.verbosity == "DEBUG" or args.raise_exception:
                raise e
            sys.exit(1)


def analysis(args):

    # Input and output filename
    infile = args.input_file
    if args.output_file is None:
        if args.output_format is None:
            raise ValueError("Extension of the output format unknown."
                             " You must either provide an output file name (with"
                             " extension) or provide it with the --output-format"
                             " argument")
        else:
            try:
                outext = extensions[args.output_format][0].lstrip(".")
            except KeyError:
                raise ValueError("No extension found for the format {}".format(args.output_format))
            outfile = infile.rsplit(".", 1)[0] + "." + outext
    else:
        outfile = args.output_file

    # Call a generic wrapper of all available conversion
    conv = Bioconvert(infile, outfile, in_fmt=args.input_format, out_fmt=args.output_format,
                      force=args.force)

    # # Users may provide information about the input file.
    # # Indeed, the input may be a FastQ file but with an extension
    # # that is not standard. For instance fq instead of fastq
    # # If so, we can use the --input-format fastq to overwrite the
    # # provided filename extension

    # no need to do this
    # if args.input_format:
    #     inext = args.input_format
    #     if not conv.inext.startswith("."):
    #         conv.inext = "." + inext

    if not conv.in_fmt:
        raise RuntimeError("convert infer the format from the extension name."
                           " So add extension to the input file name or use"
                           " --input-format option.")

    if not conv.out_fmt:
        raise RuntimeError("convert infer the format from the extension name."
                           " So add extension to the output file name or use"
                           " --output-format option.")

    # do we want to know the available methods ? If so, print info and quite
    if args.show_methods:
        print(conv.converter.available_methods)
        print("Please see http://bioconvert.readthedocs.io/en/master/"
              "references.html#bioconvert.{}.{} for details ".format(conv.name.lower(), conv.name))
        sys.exit(0)

    bioconvert.logger.info("Converting from %s to %s" % (conv.in_fmt, conv.out_fmt))

    params = {"threads": args.threads}


    if args.benchmark:
        conv.boxplot_benchmark(N=args.benchmark_N)
        import pylab
        pylab.savefig("benchmark_{}.png".format(conv.name))
    else:
        params["method"] = args.method
        conv(**params)

if __name__ == "__main__":
    main()



