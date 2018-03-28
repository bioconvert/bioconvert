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
import argparse
import json
import sys

import bioconvert
from bioconvert import ConvBase
from bioconvert.core.base import ConvMeta
from bioconvert.core.converter import Bioconvert
from bioconvert.core.decorators import get_known_dependencies_with_availability
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

        mapper = Registry()
        print("Available mapping:")
        print("==================")
        for k in sorted(mapper.get_conversions()):
            print("{} -> {}".format(k[0], k[1]))
        sys.exit(0)


class GetKnownDependenciesAction(argparse.Action):

    def __init__(self,
                 option_strings,
                 dest=argparse.SUPPRESS,
                 default=argparse.SUPPRESS,
                 help=None):
        super(GetKnownDependenciesAction, self).__init__(option_strings=option_strings,
                                                         dest=dest,
                                                         default=default,
                                                         nargs=0,
                                                         help=help)

    def __call__(self, parser, namespace, values, option_string=None):
        Registry()
        print(json.dumps(get_known_dependencies_with_availability(as_dict=True), sort_keys=True, indent=4))
        sys.exit(0)


def main(args=None):
    if args is None:
        args = sys.argv[1:]
    # Set the default level
    bioconvert.logger.level = "ERROR"

    # Changing the log level before argparse is run
    try:
        bioconvert.logger.level = args[args.index("-l") + 1]
    except:
        pass
    try:
        bioconvert.logger.level = args[args.index("--level") + 1]
    except:
        pass

    try:
        bioconvert.logger.level = args[args.index("-v") + 1]
    except:
        pass
    try:
        bioconvert.logger.level = args[args.index("--verbosity") + 1]
    except:
        pass

    allow_indirect_conversion = False
    try:
        args.index("--allow-indirect-conversion")
        allow_indirect_conversion = True
    except:
        pass
    try:
        args.index("-a")
        allow_indirect_conversion = True
    except:
        pass

    if "--version" in args:
        print("Bioconvert version {}".format(bioconvert.version))
        sys.exit(0)

    from easydev.console import purple
    if "-v" in args or "--verbosity" in args:
        print(purple("Welcome to bioconvert (bioconvert.readthedocs.io)"))

    arg_parser = argparse.ArgumentParser(prog="bioconvert",
                                         epilog=" ----    ",
                                         description="""Convertor infer the
                                         formats from the first command. We do
                                         not scan the input file. Therefore
                                         users must ensure that their input
                                         format files are properly
                                         formatted.""",
                                         usage="""
    # convert fastq to fasta
    bioconvert fastq2fasta test.fastq test.fasta
    bioconvert fastq2fasta test.fastQ test.fasta
    # if input extension is not standard, it is not a problem, format if obtained from command fastq2fasta
    bioconvert fastq2fasta test_in_fastQ.txt test.fasta

    # You may have several inputs, in which case wildcards are possible
    # Note, however, the quotes that are required
    bioconvert fastq2fasta "test*.fastq"

    # batch is also possible.
    bioconvert -m fastq2fasta "test*.fastq"

    Note the difference between the two previous commands !!


    For more information, please type:

        bioconvert --help

...
""")
    registry = Registry()
    subparsers = arg_parser.add_subparsers(help='sub-command help', dest='command', )
    max_converter_width = 2 + max([len(in_fmt) for in_fmt, _, _, _ in registry.iter_converters()])
    for in_fmt, out_fmt, converter, path in registry.iter_converters(allow_indirect_conversion):
        sub_parser_name = "{}2{}".format(in_fmt.lower(), out_fmt.lower())
        # methods = converter.available_methods if converter else []
        help_details = ""
        if converter:
            if len(converter.available_methods) <= 1:
                help_details = ""
            else:
                help_details = " (%i methods)" % len(converter.available_methods)
        elif path:
            if len(path) == 3:
                help_details = " (w/ 1 intermediate)"
            else:
                help_details = " (w/ %i intermediates)" % (len(path) - 2)
        help_text = '%sto-> %s%s' % (
            (in_fmt + ' ').ljust(max_converter_width, '-'),
            out_fmt,
            help_details,
        )
        sub_parser = subparsers.add_parser(
            sub_parser_name,
            help=help_text,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            # aliases=["{}_to_{}".format(in_fmt.lower(), out_fmt.lower()), ],
        )

        if converter:
            converter.add_argument_to_parser(sub_parser=sub_parser)
        elif path:
            for a in ConvBase.get_common_arguments():
                a.add_to_sub_parser(sub_parser)

    arg_parser.add_argument("-v", "--verbosity",
                            default=bioconvert.logger.level,
                            help="Set the outpout verbosity. Should be one of DEBUG, INFO, WARNING, ERROR, CRITICAL")
    arg_parser.add_argument("--dependency-report",
                            action=GetKnownDependenciesAction,
                            default=False,
                            help="Output all dependencies in json and exit")

    args = arg_parser.parse_args(args)

    if args.command is None:
        msg = 'No converter specified. You can list converter by doing bioconvert --help'
        arg_parser.error(msg)

    if not (getattr(args, "show_methods", False) or args.input_file):
        arg_parser.error('Either specify an input_file (<INPUT_FILE>) or ask for available methods (--show-method)')

    if not args.allow_indirect_conversion and ConvMeta.split_converter_to_extensions(args.command) not in registry:
        arg_parser.error('The conversion %s is not available directly, you have to accept that we chain converter to do'
                         ' so (--allow-indirect-conversion or -a)' % args.command)

    args.raise_exception = args.raise_exception or args.verbosity == "DEBUG"

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
            if args.raise_exception:
                raise e
            else:
                bioconvert.logger.error(e)
            sys.exit(1)


def analysis(args):
    print(vars(args))
    in_fmt, out_fmt = ConvMeta.split_converter_to_extensions(args.command)

    # do we want to know the available methods ? If so, print info and quite
    if getattr(args, "show_methods", False):
        class_converter = Registry()[(in_fmt, out_fmt)]
        print(class_converter.available_methods)
        print("Please see http://bioconvert.readthedocs.io/en/master/"
              "references.html#{} for details ".format(str(class_converter).split("'")[1]))
        if args.raise_exception:
            return
        sys.exit(0)

    # Input and output filename
    infile = args.input_file
    if args.output_file is None and infile:
        outext = ConvMeta.split_converter_to_extensions(args.command)
        outfile = infile.rsplit(".", 1)[0] + "." + outext[1].lower()
    else:
        outfile = args.output_file

    # Call a generic wrapper of all available conversion
    conv = Bioconvert(
        infile,
        outfile,
        in_fmt=in_fmt,
        out_fmt=out_fmt,
        force=args.force,
    )

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

    bioconvert.logger.info("Converting from %s to %s" % (conv.in_fmt, conv.out_fmt))

    # params = {"threads": args.threads}

    if args.benchmark:
        conv.boxplot_benchmark(N=args.benchmark_N)
        import pylab
        pylab.savefig("benchmark_{}.png".format(conv.name))
    else:
        # params["method"] = args.method
        conv(**vars(args))


if __name__ == "__main__":
    main()
