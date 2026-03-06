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
import glob
import json
import sys
import textwrap
from pathlib import Path

import colorlog
import rich_click as click

from bioconvert import logger, version
from bioconvert import ConvBase
from bioconvert.core import graph
from bioconvert.core import utils
from bioconvert.core.base import ConvMeta, ConvArg
from bioconvert.core.converter import Bioconvert
from bioconvert.core.decorators import get_known_dependencies_with_availability
from bioconvert.core.registry import Registry

_log = colorlog.getLogger(__name__)

# ── rich_click configuration ──────────────────────────────────────────────────
click.rich_click.USE_RICH_MARKUP = True
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.USE_MARKDOWN = False


# ── ConvArg → click.Parameter conversion ─────────────────────────────────────

def _conv_arg_to_click_param(conv_arg):
    """Convert a :class:`~bioconvert.core.base.ConvArg` to a :class:`click.Parameter`.

    Positional names (no leading ``-``) become :class:`click.Argument`;
    option names (leading ``-``) become :class:`click.Option`.
    """
    names = conv_arg.args_for_sub_parser
    kwargs = dict(conv_arg.kwargs_for_sub_parser)

    # Remove keys that are argparse-specific and not used by click
    kwargs.pop("output_argument", None)

    is_positional = not names[0].startswith("-")
    action = kwargs.pop("action", None)
    type_ = kwargs.pop("type", None)
    choices = kwargs.pop("choices", None)
    nargs = kwargs.pop("nargs", None)
    default = kwargs.pop("default", None)
    help_text = kwargs.pop("help", "") or ""

    # Map argparse types/choices to click equivalents
    if type_ == ConvArg.file:
        type_ = click.Path()
    elif choices:
        type_ = click.Choice(choices)

    if is_positional:
        required = nargs != "?"
        return click.Argument(
            [names[0]],
            required=required,
            default=default,
            type=type_,
        )

    # Named option
    if action == "store_true":
        return click.Option(names, is_flag=True, default=False, help=help_text)

    option_kwargs = {"help": help_text, "default": default}
    if type_:
        option_kwargs["type"] = type_
    if nargs == "+":
        option_kwargs["multiple"] = True
        # click requires a tuple/list as default when multiple=True
        if default is not None and not isinstance(default, (list, tuple)):
            option_kwargs["default"] = (default,)
        else:
            option_kwargs["default"] = () if default is None else tuple(default)
    return click.Option(names, **option_kwargs)


# ── Dynamic converter command builder ─────────────────────────────────────────

def _make_converter_command(cmd_name, converter_cls, path):
    """Build a :class:`click.RichCommand` for *cmd_name*.

    When *converter_cls* is provided the command includes all IO arguments,
    common arguments, method selection and any converter-specific additional
    arguments.  When only *path* is given (indirect conversion), only the
    basic IO and common arguments are added.
    """
    import itertools

    params = []
    seen = set()

    if converter_cls:
        arg_sources = itertools.chain(
            ConvBase.get_IO_arguments(),
            ConvBase.get_common_arguments_for_converter.__func__(converter_cls),
            converter_cls.get_additional_arguments(),
        )
    else:
        arg_sources = itertools.chain(
            ConvBase.get_IO_arguments(),
            ConvBase.get_common_arguments(),
        )

    for conv_arg in arg_sources:
        try:
            param = _conv_arg_to_click_param(conv_arg)
            # Use a stable key to de-duplicate (e.g. ('-v', '--verbosity'))
            key = tuple(conv_arg.args_for_sub_parser)
            if key not in seen:
                params.append(param)
                seen.add(key)
        except Exception as exc:
            _log.debug("Could not convert arg %s: %s", conv_arg.args_for_sub_parser, exc)

    if converter_cls:
        in_fmt, out_fmt = ConvMeta.split_converter_to_format(cmd_name)
        in_str = "/".join(in_fmt).upper()
        out_str = "/".join(out_fmt).upper()
        desc = "Convert {} to {} format. See bioconvert.readthedocs.io for details".format(
            in_str, out_str
        )
    else:
        desc = "Indirect conversion: {}".format(cmd_name)

    epilog = (
        "Bioconvert is an open source collaborative project.\n"
        "Please feel free to join us at https://github.com/bioconvert/bioconvert"
    )

    def callback(**kwargs):
        _run_analysis(cmd_name, converter_cls, path, kwargs)

    return click.RichCommand(
        name=cmd_name,
        callback=callback,
        params=params,
        help=desc,
        epilog=epilog,
    )


# ── Conversion execution ──────────────────────────────────────────────────────

def _run_analysis(converter_name, converter_cls, path, kwargs):
    """Orchestrate conversion for *converter_name* using parsed *kwargs*."""
    verbosity = kwargs.get("verbosity", "WARNING")
    logger.level = verbosity
    raise_exception = kwargs.get("raise_exception", False) or verbosity == "DEBUG"

    show_methods = kwargs.get("show_methods", False)
    input_file = kwargs.get("input_file")

    # --show-methods: print available methods and exit
    if show_methods:
        in_fmt, out_fmt = ConvMeta.split_converter_to_format(converter_name)
        class_converter = Registry()[(in_fmt, out_fmt)]
        click.echo(
            "\nMethods available for this converter ({}) are: {}".format(
                converter_name, class_converter.available_methods
            )
        )
        click.echo(
            "\nPlease see http://bioconvert.readthedocs.io/en/main/"
            "references.html#{} for details".format(str(class_converter).split("'")[1])
        )
        if not raise_exception:
            sys.exit(0)
        return

    # Require an input file unless --show-methods was given
    if input_file is None:
        raise click.UsageError(
            "Either specify an input_file (<INPUT_FILE>) or "
            "ask for available methods (--show-methods)"
        )

    # Require explicit opt-in for indirect conversions
    allow_indirect = kwargs.get("allow_indirect_conversion", False)
    if not allow_indirect and converter_cls is None:
        raise click.UsageError(
            "The conversion {} is not available directly, "
            "you have to accept that we chain converter to do so "
            "(--allow-indirect-conversion or -a)".format(converter_name)
        )

    # Handle glob / batch patterns
    if "*" in input_file or "?" in input_file:
        filenames = glob.glob(input_file)
    else:
        filenames = [input_file]

    N = len(filenames)
    for i, filename in enumerate(filenames):
        if N > 1:
            _log.info("Converting %s (%d/%d)", filename, i + 1, N)
        kwargs["input_file"] = filename
        try:
            _do_analysis(converter_name, **kwargs)
        except Exception as exc:
            if raise_exception:
                raise
            logger.error(exc)
            sys.exit(1)


def _do_analysis(converter_name, **kwargs):
    """Perform one file conversion for *converter_name*."""
    infile = kwargs["input_file"]
    outfile = kwargs.get("output_file")
    force = kwargs.get("force", False)
    threads = kwargs.get("threads")
    extra_arguments = kwargs.get("extra_arguments", "")

    # Validate tuple inputs (multi-file converters)
    if isinstance(infile, tuple):
        for file in infile:
            if not os.path.exists(file):
                if "plink" not in converter_name:
                    _log.error("Input file %s does not exist (analysis)", file)
                    sys.exit(1)

    # Auto-derive output filename when not provided
    if outfile is None and infile:
        outext = ConvMeta.split_converter_to_format(converter_name)
        if infile.split(".")[-1] in ["gz", "dsrc", "bz2"]:
            outfile = infile.rsplit(".", 1)[0].rsplit(".", 1)[0] + "." + outext[1][0].lower()
        else:
            outfile = infile.rsplit(".", 1)[0] + "." + outext[1][0].lower()

    bioconv = Bioconvert(infile, outfile, force=force, threads=threads, extra=extra_arguments)

    # Forward converter-specific extra kwargs to converter.others
    _excluded = {
        "level", "verbosity", "dependency_report", "version", "conversion_graph",
        "input_file", "output_file", "force", "extra_arguments",
        "allow_indirect_conversion", "method", "show_methods", "converter",
        "raise_exception", "batch", "benchmark", "benchmark_n", "benchmark_mode",
        "benchmark_methods", "threads", "benchmark_tag", "benchmark_save_image",
    }
    for k, v in kwargs.items():
        if k not in _excluded:
            bioconv.converter.others[k] = v

    benchmark = kwargs.get("benchmark", False)
    if benchmark:
        # benchmark_methods is a tuple when click multiple=True
        bm_methods = kwargs.get("benchmark_methods", ())
        if not bm_methods:
            bm_methods = "all"
        elif len(bm_methods) == 1 and bm_methods[0] == "all":
            bm_methods = "all"
        bioconv.converter.compute_benchmark(
            N=kwargs.get("benchmark_n", 5),
            to_include=bm_methods,
        )
        import pylab

        results = bioconv.converter.boxplot_benchmark(mode=kwargs.get("benchmark_mode", "time"))
        benchmark_tag = kwargs.get("benchmark_tag", "bioconvert")
        json_file = Path("{}.json".format(benchmark_tag))
        json_file.parent.mkdir(parents=True, exist_ok=True)

        if kwargs.get("benchmark_save_image", False):
            pylab.savefig("{}.png".format(benchmark_tag), dpi=200)
            logger.info("File %s.png created", benchmark_tag)

        with open(json_file, "w") as fout:
            json.dump(results, fout, indent=True, sort_keys=True)
            logger.info("Saved results in %s", json_file)
    else:
        bioconv(**kwargs)


# ── Custom Click Group with dynamic subcommand loading ────────────────────────

class BioconvertCLI(click.RichGroup):
    """Click group that dynamically builds subcommands from the converter registry.

    Supports **implicit mode**: if the first argument looks like a filename
    (contains a ``.``) the converter name is automatically inferred from the
    file extensions before dispatching.
    """

    def __init__(self, **attrs):
        super().__init__(**attrs)
        self._registry = None

    @property
    def registry(self):
        if self._registry is None:
            logger.level = "ERROR"
            self._registry = Registry()
            logger.level = "WARNING"
        return self._registry

    # ------------------------------------------------------------------
    # Implicit-mode pre-processing
    # ------------------------------------------------------------------

    def parse_args(self, ctx, args):
        """Intercept args before Click processes them.

        When the first positional argument looks like a filename (contains
        a ``.``) and is not a known converter name the converter is
        auto-detected from file extensions (*implicit mode*).

        Also pre-detects ``-a / --allow-indirect-conversion`` so that
        :meth:`get_command` and :meth:`list_commands` can use it before the
        group's own option parsing runs.
        """
        args = list(args)

        # Pre-detect allow-indirect-conversion so get_command can use it
        # even when the flag appears after the subcommand name.
        allow_indirect = "-a" in args or "--allow-indirect-conversion" in args
        obj = ctx.ensure_object(dict)
        obj["allow_indirect_conversion"] = allow_indirect

        if args and not args[0].startswith("-"):
            reg = self.registry
            if args[0].lower() not in list(reg.get_converters_names()) and "." in args[0]:
                args = self._handle_implicit_mode(args, reg)
        return super().parse_args(ctx, args)

    def _handle_implicit_mode(self, args, registry):
        """Detect converter from file extensions and inject it as the first argument."""
        if not os.path.exists(args[0]):
            _log.error("First input file %s does not exist", args[0])
            sys.exit(1)

        # Collect positional arguments (stop at the first option flag)
        positional = []
        for arg in args:
            if arg.startswith("-"):
                break
            positional.append(arg)

        filenames = positional
        exts = [utils.get_extension(x, remove_compression=True) for x in filenames]

        # Try all possible input/output splits
        L = len(exts)
        converter_list = []
        in_ext = out_ext = None
        for i in range(1, L):
            in_ext = tuple(exts[0:i])
            out_ext = tuple(exts[i:])
            try:
                converter_list.extend(registry.get_ext((in_ext, out_ext)))
            except KeyError:
                pass

        # Special case: same extension with different compression
        if not converter_list and len(exts) >= 2 and exts[0] == exts[1]:
            exts_with_comp = [utils.get_extension(x, remove_compression=False) for x in filenames]
            in_ext_c, out_ext_c = exts_with_comp[0], exts_with_comp[1]
            comps = ["gz", "dsrc", "bz2"]
            if in_ext_c in comps and out_ext_c in comps:
                try:
                    converter_list.extend(registry.get_ext(((in_ext_c,), (out_ext_c,))))
                except KeyError:
                    pass

        if not converter_list:
            in_label = exts[0] if exts else "?"
            out_label = exts[-1] if len(exts) > 1 else "?"
            msg = "\nBioconvert does not support conversion {} -> {}. \n\n".format(
                in_label, out_label
            )
            # Suggest indirect path if available
            try:
                _path = registry._path_dict_ext[in_label][out_label]
                in_name, int_name, out_name = _path
                a = registry._ext_registry[in_name, int_name][0].__name__.split("2")[0]
                b = registry._ext_registry[int_name, out_name][0].__name__.split("2")[1]
                convname = "2".join([a, b]).lower()
                msg += "\n".join(
                    textwrap.wrap(
                        "Note, however, that an indirect conversion through "
                        "an intermediate format is possible. To do so, use the "
                        "-a option and be explicit about the conversion type."
                    )
                )
                msg += "\n\n    bioconvert --help -a\n\n"
                msg += "The command should probably be:\n\n"
                msg += "    bioconvert {} {} -a\n\n".format(convname, " ".join(positional))
            except (KeyError, ValueError):
                pass
            _log.error(msg)
            sys.exit(1)
        elif len(converter_list) == 1:
            converter_name = converter_list[0].__name__.lower()
        else:
            _log.error(
                "Ambiguous extension.\nYou must specify the right conversion. "
                "Please choose from:\n\n%s",
                "\n".join([c.__name__.lower() for c in converter_list]),
            )
            sys.exit(1)

        return [converter_name] + args

    # ------------------------------------------------------------------
    # Dynamic subcommand discovery
    # ------------------------------------------------------------------

    def list_commands(self, ctx):
        """Return sorted list of available converter names."""
        allow_indirect = ctx.params.get("allow_indirect_conversion", False)
        if not allow_indirect and ctx.obj:
            allow_indirect = ctx.obj.get("allow_indirect_conversion", False)
        commands = []
        for in_fmt, out_fmt, _conv, _path in self.registry.iter_converters(allow_indirect):
            in_fmt = ConvBase.lower_tuple(in_fmt)
            out_fmt = ConvBase.lower_tuple(out_fmt)
            commands.append("{}2{}".format("_".join(in_fmt), "_".join(out_fmt)))
        return sorted(commands)

    def get_command(self, ctx, cmd_name):
        """Dynamically build and return the click Command for *cmd_name*."""
        cmd_name = cmd_name.lower()

        try:
            in_fmt, out_fmt = ConvMeta.split_converter_to_format(cmd_name)
        except TypeError:
            return None

        allow_indirect = ctx.params.get("allow_indirect_conversion", False)
        if not allow_indirect and ctx.obj:
            allow_indirect = ctx.obj.get("allow_indirect_conversion", False)
        registry = self.registry

        converter_cls = None
        path = None

        if (in_fmt, out_fmt) in registry:
            converter_cls = registry[(in_fmt, out_fmt)]
        elif allow_indirect:
            try:
                path = registry.conversion_path(in_fmt, out_fmt)
            except Exception:
                pass
            if not path:
                return None
        else:
            return None

        return _make_converter_command(cmd_name, converter_cls, path)


# ── Main CLI group ─────────────────────────────────────────────────────────────

_EPILOG = """
Bioconvert contains tens of converters whose list is available as follows:

    bioconvert --help

Each conversion has its own sub-command and dedicated help. For instance:

    bioconvert fastq2fasta --help

Because the subcommand contains the format, extensions are not important
for the conversion itself. This would convert the test.txt file (fastq
format) into a fasta file:

    bioconvert fastq2fasta test.txt test.fasta

If you use known extensions, the converter may be omitted:

    bioconvert test.fastq test.fasta

Users must ensure that their input format files are properly formatted.

If there is a conversion from A to B and another for B to C, you can also
perform indirect conversion using -a argument (experimental). This command
shows all possible indirect conversions:

    bioconvert --help -a

Please visit http://bioconvert.readthedocs.org for more information.
Would you wish to help, please join our open source collaborative project
at https://github.com/bioconvert/bioconvert
"""


@click.group(
    cls=BioconvertCLI,
    invoke_without_command=True,
    epilog=_EPILOG,
    context_settings={"help_option_names": ["--help"]},
)
@click.option(
    "-v", "--verbosity",
    default=logger.level,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
    help="Set the output verbosity.",
)
@click.option(
    "-l", "--level",
    default=logger.level,
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]),
    help="Set the output verbosity. Same as --verbosity.",
    hidden=True,
)
@click.option(
    "--dependency-report",
    "dependency_report",
    is_flag=True,
    default=False,
    help="Output all bioconvert dependencies in JSON and exit.",
)
@click.option(
    "-a", "--allow-indirect-conversion",
    "allow_indirect_conversion",
    is_flag=True,
    default=False,
    help="Show / allow indirect conversions (labelled as intermediate).",
)
@click.option(
    "--version",
    "show_version",
    is_flag=True,
    default=False,
    help="Show version and exit.",
)
@click.option(
    "--conversion-graph",
    "conversion_graph",
    default=None,
    type=click.Choice(["cytoscape", "cytoscape-all"]),
    help="Output the conversion graph in the specified format.",
)
@click.pass_context
def cli(ctx, verbosity, level, dependency_report, allow_indirect_conversion,
        show_version, conversion_graph):
    """Bioconvert: convert between bioinformatics file formats.

    Use ``bioconvert CONVERTER --help`` for help on a specific converter.
    """
    logger.level = verbosity or level

    if show_version:
        click.echo("{}".format(version))
        sys.exit(0)

    if dependency_report:
        click.echo(
            json.dumps(
                get_known_dependencies_with_availability(as_dict=True),
                sort_keys=True,
                indent=4,
            )
        )
        sys.exit(0)

    if conversion_graph:
        if conversion_graph.startswith("cytoscape"):
            all_converter = conversion_graph == "cytoscape-all"
            click.echo(
                json.dumps(
                    graph.create_graph_for_cytoscape(all_converter=all_converter),
                    indent=4,
                )
            )
        sys.exit(0)

    # When invoked without a subcommand and no special flags, raise a usage error
    if ctx.invoked_subcommand is None:
        raise click.UsageError(
            "No converter specified. "
            "You can list all converters by using:\n\n\tbioconvert --help"
        )


# ── Public entry point ────────────────────────────────────────────────────────

def _find_sub_command(args):
    """Return the first positional argument that could be a converter name."""
    skip_next = False
    _option_with_value = {"-v", "--verbosity", "-l", "--level", "--conversion-graph"}
    for arg in args:
        if skip_next:
            skip_next = False
            continue
        if arg in _option_with_value:
            skip_next = True
            continue
        if not arg.startswith("-"):
            return arg
    return None


def _suggest_close_matches(sub_command, registry):
    """Return up to 5 converter names within Levenshtein distance 3 of *sub_command*."""
    from bioconvert.core.levenshtein import wf_levenshtein as lev

    conversions = []
    for in_fmt, out_fmt, _conv, _path in registry.iter_converters(False):
        in_fmt = ConvBase.lower_tuple(in_fmt)
        out_fmt = ConvBase.lower_tuple(out_fmt)
        cname = "{}2{}".format("_".join(in_fmt), "_".join(out_fmt))
        conversions.append((lev(cname, sub_command), cname))

    return sorted([x for x in conversions if x[0] <= 3])[:5]


def main(args=None):
    """Entry point for the ``bioconvert`` CLI command."""
    if args is None:
        args = sys.argv[1:]
    args = list(args)

    # Pre-set log level early so it applies during registry initialisation
    for flag in ("-v", "--verbosity", "-l", "--level"):
        try:
            idx = args.index(flag)
            if idx + 1 < len(args):
                logger.level = args[idx + 1]
        except (ValueError, IndexError):
            pass

    # No arguments at all → exit 1 (consistent with previous behaviour)
    if not args:
        _log.error("Please provide at least some arguments. See --help")
        sys.exit(1)

    try:
        # standalone_mode=False: click returns the exit code for clean exits
        # (e.g. --help triggers Exit(0) → returns 0) and raises ClickException
        # for usage errors.
        result = cli(args=args, standalone_mode=False)
    except SystemExit:
        # Explicit sys.exit() calls (e.g. from --version, --show-methods)
        # propagate unchanged.
        raise
    except click.Abort:
        sys.exit(1)
    except click.UsageError as exc:
        # A usage error might be a typo in the converter name; attempt to
        # suggest close matches before falling back to exit 2.
        logger.level = "ERROR"
        registry = Registry()
        logger.level = "WARNING"

        sub_command = _find_sub_command(args)

        if sub_command is not None:
            matches = _suggest_close_matches(sub_command, registry)

            # Sub-command exists but something else is wrong → keep exit 2
            if matches and sub_command in [x[1] for x in matches]:
                _log.error(str(exc))
                sys.exit(2)

            if matches:
                msg = (
                    "\n\nThe converter {}() was not found. It may be a typo.\n"
                    "Here is a list of possible matches: {} ...\n"
                    "You may also add the -a argument to enforce a transitive "
                    "conversion. The whole list is available using\n\n"
                    "    bioconvert --help -a\n"
                )
                _log.error(msg.format(sub_command, ", ".join([v for _, v in matches])))
                sys.exit(1)

        _log.error(str(exc))
        sys.exit(2)
    else:
        # With standalone_mode=False, click returns the exit code when a clean
        # exit is requested (e.g. --help triggers Exit(0) → cli() returns 0).
        if result is not None:
            sys.exit(result)


if __name__ == "__main__":
    main()
