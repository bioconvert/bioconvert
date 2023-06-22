###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
""".. rubric:: Standalone application dedicated to conversion"""
import os
import argparse
import glob
import json
import sys
import colorlog
import textwrap
import json
from pathlib import Path

from bioconvert import logger, version
from bioconvert import ConvBase
from bioconvert.core import graph
from bioconvert.core import utils
from bioconvert.core.base import ConvMeta
from bioconvert.core.converter import Bioconvert
from bioconvert.core.decorators import get_known_dependencies_with_availability
from bioconvert.core.registry import Registry

_log = colorlog.getLogger(__name__)


def error(msg):
    _log.error(msg)
    sys.exit(1)


class ParserHelper:
    """
    convenient variable to check implicit/explicit mode and
    get information about the arguments
    """

    def __init__(self, args):
        # We do not need boring lengthy warning at that stage.
        logger.level = "ERROR"
        registry = Registry()
        logger.level = "WARNING"

        self.args = args[:]
        if len(self.args) == 0:
            error("Please provide at least some arguments. See --help")

        if args[0].lower() not in list(registry.get_converters_names()) and "." in args[0]:
            self.mode = "implicit"
            # we should have at least 2 entries in implicit mode and the first
            # input filename must exists (1 to 1 or many to 1)
            if len(args) < 2:
                error("In implicit mode, you must define your input and output file (only 1 provided)")
        else:
            self.mode = "explicit"
        _log.debug("parsing mode {}".format(self.mode))

    def get_filelist(self):
        """return list of the input files"""
        positional = []
        for arg in self.args:
            if arg.startswith("-"):
                break
            else:
                positional.append(arg)

        # remaining arguments are files
        if self.mode == "implicit":
            return positional
        else:
            return positional[1:]


def main(args=None):


    # used later on
    registry = Registry()

    if args is None:
        args = sys.argv[1:]

    # convenient variable to check implicit/explicit mode and
    # get information about the arguments.
    ph = ParserHelper(args)

    if not len(sys.argv) == 1:

        if ph.mode == "implicit":

            # Check that the input file exists
            # Fixes https://github.com/bioconvert/bioconvert/issues/204
            if os.path.exists(args[0]) is False:
                _log.error("First input file {} does not exist".format(args[0]))
                sys.exit(1)

            # list of filenames from which we get the extensions
            filenames = ph.get_filelist()
            exts = [utils.get_extension(x, remove_compression=True) for x in filenames]

            # We need to get the corresponding converter if any.

            # We assume that the input formats are ordered alphabetically
            # (bioconvert API).
            # For instance fasta,qual to fastq can be
            # found but qual,fasta to fastq cannot. Indeed, in more complex
            # cases such as a,b -> c,d we cannot know whether there are 1 or 3
            # inputs. This would require extra code here below
            try:
                L = len(exts)
                converter = []
                # if input is a,b,c,d we want to try a->(b,c,d) and
                # (a,b)->(c,d) and (a,b,c)-> c so L-1 case
                for i in range(1, L):
                    in_ext = tuple(exts[0:i])
                    out_ext = tuple(exts[i:])
                    try:
                        converter.extend(registry.get_ext((in_ext, out_ext)))
                    except KeyError:
                        pass
            except KeyError:
                converter = []

            # For 1-to-1, if the extensions are identical but different
            # compression, this means we just want to decompress and
            # re-compress in another format.
            if not converter and (exts[0] == exts[1]):
                exts_with_comp = [utils.get_extension(x, remove_compression=False) for x in filenames]
                in_ext, out_ext = exts_with_comp[0], exts_with_comp[1]
                comps = ["gz", "dsrc", "bz2"]
                if in_ext in comps and out_ext in comps:
                    converter.extend(registry.get_ext(((in_ext,), (out_ext,))))

            # if no converter is found, print information
            if not converter:
                msg = "\nBioconvert does not support conversion {} -> {}. \n\n"
                msg = msg.format(in_ext, out_ext)

                # maybe it is an indirect conversion ? let us look at the
                # digraph
                try:
                    _path = registry._path_dict_ext[in_ext][out_ext]
                    # Here, we have a transitive list of tuples to go from A to C
                    # example from fq to clustal returns:
                    # [('fq',), ('fa',), ('clustal',)]
                    # If we naively build the converter from those names
                    # (fq2clustal), this is a non official converter name. The
                    # official one is fastq2clustal, so we need some hack here:
                    in_name, int_name, out_name = _path
                    a = registry._ext_registry[in_name, int_name][0].__name__.split("2")[0]
                    b = registry._ext_registry[int_name, out_name][0].__name__.split("2")[1]

                    convname = "2".join([a, b]).lower()

                    msg += "\n".join(
                        textwrap.wrap(
                            "Note, however, that an indirect conversion through"
                            " an intermediate format is possible for your input and "
                            " output format. To do so, you need to use the -a option "
                            " and be explicit about the type of conversion. To get "
                            " the list of possible direct and indirect conversion, "
                            " please use:\n\n"
                        )
                    )
                    msg += "\n\n    bioconvert --help -a\n\n"
                    msg += "For help and with your input/output most probably"
                    msg += "the command should be: \n\n    bioconvert {} {} -a\n\n ".format(
                        convname, " ".join(ph.get_filelist())
                    )
                except KeyError:
                    pass  # not converter found in the path
                error(msg)

            # if the ext_pair matches a single converter
            elif len(converter) == 1:
                args.insert(0, converter[0].__name__.lower())
            # if the ext_pair matches multiple converters
            else:
                _log.error(
                    "Ambiguous extension.\n"
                    "You must specify the right conversion  Please "
                    "choose a conversion from: \n\n"
                    "{}".format("\n".join([c.__name__.lower() for c in converter]))
                )
                sys.exit(1)

    # Set the default level
    logger.level = "ERROR"

    # Changing the log level before argparse is run
    try:
        logger.level = args[args.index("-l") + 1]
    except:
        pass
    try:
        logger.level = args[args.index("--level") + 1]
    except:
        pass
    try:
        logger.level = args[args.index("-v") + 1]
    except:
        pass
    try:
        logger.level = args[args.index("--verbosity") + 1]
    except:
        pass


    # if there is the ability to convert from A to B to C, we must set
    # the option -a (--allow_indirect_conversion)
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

    # Now, the instanciation of the main bioconvert user interface
    arg_parser = argparse.ArgumentParser(
        prog="bioconvert",
        description="",
        # ""Convertor infer the
        # formats from the first command. We do
        # not scan the input file. Therefore
        # users must ensure that their input
        # format files are properly
        # formatted.""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Bioconvert contains tens of converters whose list is available as follows:

    bioconvert --help

Each conversion has its own sub-command and dedicated help. For instance:

    bioconvert fastq2fasta --help

Because the subcommand contains the format, extensions are not important
for the conversion itself. This would convert the test.txt file (fastq
format) into a fasta file:

    bioconvert fastq2fasta test.txt test.fasta

If you use known extensions, the converter may be omitted::

    bioconvert test.fastq test.fasta

Users must ensure that their input format files are properly formatted.

If there is a conversion from A to B and another for B to C, you can also
perform indirect conversion using -a argument (experimental). This command
shows all possible indirect conversions:

    bioconvert --help -a

Please visit http://bioconvert.readthedocs.org for more information about the
project or formats available. Would you wish to help, please join our open 
source collaborative project at https://github/bioconvert/bioconvert
""",
    )

    subparsers = arg_parser.add_subparsers(
        help="sub-command help",
        dest="converter",
    )

    max_converter_width = 2 + max([len(in_fmt) for in_fmt, _, _, _ in registry.iter_converters()])

    def sorting_tuple_string(item):
        if type(item) is tuple:
            return item[0][0]
        if type(item) is str:
            return item[0]

    # show all possible conversion including indirect conversion
    for in_fmt, out_fmt, converter, path in sorted(
        registry.iter_converters(allow_indirect_conversion), key=sorting_tuple_string
    ):
        in_fmt = ConvBase.lower_tuple(in_fmt)
        in_fmt = ["_".join(in_fmt)]

        out_fmt = ConvBase.lower_tuple(out_fmt)
        out_fmt = ["_".join(out_fmt)]

        sub_parser_name = "{}2{}".format("_".join(in_fmt), "_".join(out_fmt))

        if converter:
            link_char = "-"
            if len(converter.available_methods) < 1:
                help_details = " (no available methods please see the doc" " for install the necessary libraries) "
            else:
                help_details = " (%i methods)" % len(converter.available_methods)
        else:  # if path:
            link_char = "~"
            if len(path) == 3:
                help_details = " (w/ 1 intermediate)"
            else:
                help_details = " (w/ %i intermediates)" % (len(path) - 2)

        help_text = "{}to{}> {}{}".format(
            ("_".join(in_fmt) + " ").ljust(max_converter_width, link_char),
            link_char,
            ("_".join(out_fmt)),
            help_details,
        )
        sub_parser = subparsers.add_parser(
            sub_parser_name,
            help=help_text,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            # aliases=["{}_to_{}".format(in_fmt.lower(), out_fmt.lower()), ],
            epilog="""Bioconvert is an open source collaborative project. 
Please feel free to join us at https://github/biokit/bioconvert
""",
        )
        if converter:
            converter.add_argument_to_parser(sub_parser=sub_parser)
        elif path:
            for a in ConvBase.get_IO_arguments():
                a.add_to_sub_parser(sub_parser)
            for a in ConvBase.get_common_arguments():
                a.add_to_sub_parser(sub_parser)

    # arguments when no explicit conversion provided.

    arg_parser.add_argument(
        "-v",
        "--verbosity",
        default=logger.level,
        help="Set the outpout verbosity. Same as --level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )

    arg_parser.add_argument(
        "-l",
        "--level",
        default=logger.level,
        help="Set the outpout verbosity. Same as --verbosity",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
    )

    arg_parser.add_argument(
        "--dependency-report",
        action="store_true",
        default=False,
        help="Output all bioconvert dependencies in json and exit",
    )

    arg_parser.add_argument(
        "-a",
        "--allow-indirect-conversion",
        action="store_true",
        help="Show all possible indirect conversions " "(labelled as intermediate)",
    )

    arg_parser.add_argument("--version", action="store_true", default=False, help="Show version")

    arg_parser.add_argument(
        "--conversion-graph",
        nargs="?",
        default=None,
        choices=[
            "cytoscape",
            "cytoscape-all",
        ],
    )

    # FIXME why is this a try ?
    # Possible to check whether the subcommand is valid or not
    try:
        args = arg_parser.parse_args(args)
    except SystemExit as err:

        # parsing ask to stop, maybe a normal exit
        if err.code == 0:
            raise err

        # Parsing failed, trying to guess converter
        from bioconvert.core.levenshtein import wf_levenshtein as lev

        sub_command = None
        args_i = 0
        while sub_command is None and args_i < len(args):
            if args[args_i][0] != "-" and (
                args_i == 0
                or args[args_i - 1] != "-v"
                and args[args_i - 1] != "--verbose"
                and args[args_i - 1] != "--conversion-graph"
            ):
                sub_command = args[args_i]
            args_i += 1

        if sub_command is None:
            # No sub_command found, so letting the initial exception be risen
            raise err

        conversions = []
        for in_fmt, out_fmt, converter, path in registry.iter_converters(allow_indirect_conversion):
            in_fmt = ConvBase.lower_tuple(in_fmt)
            in_fmt = ["_".join(in_fmt)]
            out_fmt = ConvBase.lower_tuple(out_fmt)
            out_fmt = ["_".join(out_fmt)]
            conversion_name = "{}2{}".format("_".join(in_fmt), "_".join(out_fmt))
            conversions.append((lev(conversion_name, sub_command), conversion_name))

        conversions = [x for x in conversions if x[0] <= 3]
        matches = sorted(conversions)[:5]

        # It could be a typo elsewhere, a real error in the command line such as
        # --verbosity INFFO. If so, the subcommand may be correct and therefore
        # we should stop here to not confuse the users
        if len(matches) > 0 and sub_command in [x[1] for x in matches]:
            raise err

        if len(matches) == 0:
            # sub_command was ok, problem comes from elsewhere
            msg = "sub command {} not found and no close candidate found. "
            msg = msg.format(sub_command)
        else:
            msg = (
                "\n\nThe converter {}() was not found. It may be a typo.\n"
                + "Here is a list of possible matches: {} ... "
                + "\nYou may also add the -a argument to enfore a "
                + "transitive conversion. The whole list is available using\n\n"
                + "    bioconvert --help -a \n"
            )
            msg = msg.format(sub_command, ", ".join([v for _, v in matches]))
        _log.error(msg)
        sys.exit(1)

    if args.version:
        print("{}".format(version))
        sys.exit(0)

    if args.dependency_report:
        print(
            json.dumps(
                get_known_dependencies_with_availability(as_dict=True),
                sort_keys=True,
                indent=4,
            )
        )
        sys.exit(0)

    if args.conversion_graph:
        if args.conversion_graph.startswith("cytoscape"):
            all_converter = args.conversion_graph == "cytoscape-all"
            print(
                json.dumps(
                    graph.create_graph_for_cytoscape(all_converter=all_converter),
                    indent=4,
                )
            )
        sys.exit(0)

    if args.converter is None:
        msg = "No converter specified. "
        msg += "You can list all converters by using:\n\n\tbioconvert --help"
        arg_parser.error(msg)

    if getattr(args, "show_methods", False) is False and args.input_file is None:
        arg_parser.error("Either specify an input_file (<INPUT_FILE>) or " "ask for available methods (--show-methods)")

    # do we want to know the available methods ? If so, print info and quit
    if getattr(args, "show_methods", False) is True:
        in_fmt, out_fmt = ConvMeta.split_converter_to_format(args.converter)
        class_converter = Registry()[(in_fmt, out_fmt)]
        print(
            "\nMethods available for this converter ({}) are: {}".format(
                args.converter, class_converter.available_methods
            )
        )
        print(
            "\nPlease see http://bioconvert.readthedocs.io/en/main/"
            "references.html#{} for details ".format(str(class_converter).split("'")[1])
        )
        if args.raise_exception:
            return
        sys.exit(0)

    if not args.allow_indirect_conversion and ConvMeta.split_converter_to_format(args.converter) not in registry:

        arg_parser.error(
            "The conversion {} is not available directly, "
            "you have to accept that we chain converter to do"
            " so (--allow-indirect-conversion or -a)".format(args.converter)
        )

    args.raise_exception = args.raise_exception or args.verbosity == "DEBUG"

    # Set the logging level
    logger.level = args.verbosity
    # Figure out whether we have several input files or not
    if "*" in args.input_file or "?" in args.input_file:
        filenames = glob.glob(args.input_file)
    else:
        filenames = [args.input_file]

    N = len(filenames)
    for i, filename in enumerate(filenames):
        if N > 1:
            _log.info("Converting {} ({}/{})".format(filename, i + 1, N))
        args.input_file = filename
        try:
            analysis(args)
        except Exception as e:
            if args.raise_exception:
                raise e
            else:
                logger.error(e)
            sys.exit(1)


def analysis(args):
    in_fmt, out_fmt = ConvMeta.split_converter_to_format(args.converter)

    # Input and output filename
    infile = args.input_file

    # Check that the input file exists
    # Fixes https://github.com/bioconvert/bioconvert/issues/204
    if type(infile) is tuple:
        for file in infile:
            if os.path.exists(file) is False:

                # Some convertors uses prefix instead of filename. We could have
                # ambiguities: if we use a prefix without extension,
                # we could be confused with the convertor name. This is true
                # for the plink families
                if "plink" in args.converter:
                    pass
                else:
                    _log.error("Input file {} does not exist (analysis)".format(file))
                    sys.exit(1)

    if args.output_file is None and infile:
        outext = ConvMeta.split_converter_to_format(args.converter)
        if infile.split(".")[-1] in ["gz", "dsrc", "bz2"]:
            # get rid of extension gz/dsrc/bz2
            outfile = infile.rsplit(".", 1)[0]
            # get rid of extension itself
            outfile = outfile.rsplit(".", 1)[0]
            outfile += "." + outext[1][0].lower()

        else:
            outfile = infile.rsplit(".", 1)[0] + "." + outext[1][0].lower()
    else:
        outfile = args.output_file

    # check whether a valid --thread option was provided
    if "threads" in args:
        threads = args.threads
    else:
        threads = None

    # default will be ""
    if "extra_arguments" in args:
        extra_arguments = args.extra_arguments

    # Call a generic wrapper of all available conversion
    bioconv = Bioconvert(infile, outfile, force=args.force, threads=threads, extra=extra_arguments)

    for k, v in args.__dict__.items():
        if k not in [
            "level",
            "verbosity",
            "dependency_report",
            "version",
            "conversion_graph",
            "input_file",
            "output_file",
            "force",
            "extra_arguments",
            "allow_indirect_conversion",
            "method",
            "show_methods",
            "converter",
            "raise_exception",
            "batch",
            "benchmark",
            "benchmark_N",
            "benchmark_mode",
            "benchmark_methods",
        ]:
            bioconv.converter.others[k] = v

    if args.benchmark:

        bioconv.converter.compute_benchmark(N=args.benchmark_N, to_include=args.benchmark_methods)
        import pylab

        results = bioconv.converter.boxplot_benchmark(mode=args.benchmark_mode)

        json_file = Path(f"{args.benchmark_tag}.json")
        json_file.parent.mkdir(parents=True, exist_ok=True)

        if args.benchmark_save_image:
            pylab.savefig(f"{args.benchmark_tag}.png", dpi=200)
            logger.info(f"File {args.benchmark_tag}.png created")

        with open(json_file, "w") as fout:
            json.dump(results, fout, indent=True, sort_keys=True)
            logger.info(f"Saved results in {json_file}")
    else:
        bioconv(**vars(args))


if __name__ == "__main__":
    main()
