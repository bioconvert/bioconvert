"""Convert :term:`CSV` format to :term:`TSV` file"""

import colorlog

from bioconvert import tsv2csv
from bioconvert.core.decorators import requires, requires_nothing

logger = colorlog.getLogger(__name__)


class CSV2TSV(tsv2csv.TSV2CSV):
    """Bam2Json converter

    Convert bam file to json file.
    """
    _default_method = "python"
    DEFAULT_IN_SEP = ','
    DEFAULT_OUT_SEP = '\t'
    DEFAULT_LINE_TERMINATOR = '\n'

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        """
        super().__init__(infile, outfile)

    @requires_nothing
    def _method_python(
            self,
            in_sep=DEFAULT_IN_SEP,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            *args, **kwargs):
        """
        do the conversion :term`CSV` -> :term:'TSV` using standard Python modules
        """
        super()._method_python(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)

    @requires_nothing
    def _method_python_v2(
            self,
            in_sep=DEFAULT_IN_SEP,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            *args, **kwargs):
        """
        do the conversion :term`TSV` -> :term:'CSV` using csv module to read, and writing directly into the file. Note
        that his method can't escape nor quote output char
        """
        super()._method_python_v2(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)

    @requires(python_library="pandas")
    def _method_panda(
            self,
            in_sep=DEFAULT_IN_SEP,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            *args, **kwargs):
        """
        do the conversion :term`CSV` -> :term:'TSV` using Panda modules
        """
        super()._method_panda(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)
