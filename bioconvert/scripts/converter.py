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
        verbo_nb = sum([1 for opt in sys.argv if opt.startswith('--verb')])
        verbosity = v_nb + verbo_nb

        bioconvert.logger_set_level(max(10, 30 - (10 * verbosity)))

        mapper = Registry()
        print("Available mapping:")
        print("==================")
        for k in sorted(mapper.get_conversions()):
            print("{} -> {}".format(k[0], k[1]))
        sys.exit(0)

def main(args=None):

    if args is None:
        args = sys.argv[:]

    from easydev.console import purple, underline
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

""")
    arg_parser.add_argument("input_file",
            default=None,
            help="The path to the file to convert.")
    arg_parser.add_argument("output_file", nargs="?",
            default=None,
            help="The path where the result will be stored.")

    arg_parser.add_argument("-f", "--formats",
                            action=ConvAction,
                            default=False,
                            help="Display available formats and exit.")
    arg_parser.add_argument("-v", "--verbosity",
                            action="count",
                            default=0,
                            help="Set the outpout verbosity.")
    arg_parser.add_argument("-i", "--input-format",
                            default=None,
                            help="Provide the input format. Check the --formats to see valid input name")
    arg_parser.add_argument("-o", "--output-format",
                            default=None,
                            help="Provide the output format. Check the --formats to see valid input name")
    arg_parser.add_argument("-x", "--threads",
                            default=None,
                            help="Number of threads. Depends on the underlying tool")
    arg_parser.add_argument("-m", "--batch",
                            default=False, action="store_true",
                            help="for batch effect")

    arg_parser.add_argument("-c", "--method",
                            default=None,
                            help="A converter may have several methods")

    arg_parser.add_argument("-s", "--show-methods",
                            default=False,
                            action="store_true",
                            help="A converter may have several methods")

    args = arg_parser.parse_args()

    # Set the logging level
    args.verbosity = max(10, 30 - (10 * args.verbosity))
    bioconvert.logger_set_level(args.verbosity)
    _log = colorlog.getLogger('bioconvert')


    # Figure out whether we have several input files or not
    # Are we in batch mode ? 
    import glob
    if args.batch:
        filenames = glob.glob(args.input_file)
    else:
        filenames = [args.input_file]
    print("-----------")
    print(filenames)

    for filename in filenames:
        print(filename)
        args.input_file = filename
        analysis(args)


    #_log.info("Done")


def analysis(args):
    mapper = Registry()
    _log = colorlog.getLogger('bioconvert')

    # Input and output filename
    infile = args.input_file
    if args.output_file is None:
        if args.output_format is None:
            raise ValueError("Extension of the output format unknown."
                    " You must either provide an output file name (with"
                    " extension) or provide it zith the --output-format"
                    " argument")
        else:
            outfile = infile.rsplit(".",1)[0] + "." + args.output_format
    else:
        outfile = args.output_file

    # Users may provide information about the input file.
    # Indeed, the input may be a FastQ file but with an extension
    # that is not standard. For instance fq instead of fastq
    # If so, we can use the --input-format fastq to overwrite the
    # provided filename extension
    inext = os.path.splitext(infile)[-1]
    outext = os.path.splitext(outfile)[-1]

    if args.input_format:
        inext = args.input_format
        if not inext.startswith("."):
            inext = "." + inext

    if not inext:
        raise RuntimeError("convert infer the format from the extension name."
                           " So add extension to the input file name or use"
                           " --input-format option.")

    if not outext:
        raise RuntimeError("convert infer the format from the extension name."
                           " So add extension to the output file name or use"
                           " --outut-format option.")

    # From the input parameters 1 and 2, we get the module name
    try:
        _log.info("Input: {}".format(inext))
        _log.info("Output: {}".format(outext))
        class_converter = mapper[(inext, outext)]
    except KeyError:
        print(mapper)
        print(inext)
        print(outext)

        # Is the module name available in biokit ? If not, let us tell the user
        msg = "Request input format ({}) to output format (({}) is not available in converters"
        _log.critical(msg.format(inext, outext))
        _log.critical("Use --formats to know the available formats")
        sys.exit(1)

    # If the module exists, it is part of the MapperRegitry dictionary and
    # we should be able to import it dynamically, create the class and call
    # the instance
    _log.info("Converting from {} to {}".format(inext, outext))

    # Prepare some user arguments
    params = {"threads": args.threads}
    if args.method:
        params["method"] = args.method


    # Call the class method that does the real work
    convert = class_converter(infile, outfile)

    # do we want to know the available methods ? If so, print info and quite
    if args.show_methods:
        print(convert.available_methods)
        print("Please see http://bioconvert.readthedocs.io/en/master/references.html#bioconvert.{}.{} "
              "for details ".format(class_converter.__name__.lower(),class_converter.__name__))
        sys.exit(0)

    convert(**params)



if __name__ == "__main__":
    main()



