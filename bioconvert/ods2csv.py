"""Convert :term:`XLS` format to :term:`CSV` file"""
import csv

import colorlog

from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires, requires_nothing

from bioconvert import ConvBase

logger = colorlog.getLogger(__name__)


class ODS2CSV(ConvBase):
    """Convert :term:`XLS` file into :term:`CSV` file

    """
    _default_method = "pyexcel"
    DEFAULT_OUT_SEP = ','
    DEFAULT_LINE_TERMINATOR = '\n'

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        """
        super().__init__(infile, outfile)

    @requires(python_libraries=["pyexcel", "pyexcel-ods3"])
    def _method_pyexcel(
            self,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            sheet_name=0,
            *args, **kwargs):
        """
        do the conversion :term`XLS` -> :term:'CSV` using Panda modules

        """
        import pyexcel
        with open(self.outfile, "w") as out_stream:
            writer = csv.writer(out_stream, delimiter=out_sep, lineterminator=line_terminator)
            first_row = True
            for row in pyexcel.get_records(file_name=self.infile):
                if first_row:
                    writer.writerow([k for k, v in row.items()])
                    first_row = False
                writer.writerow([v for k, v in row.items()])

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names=["--sheet-name", ],
            default=0,
            help="The name or id of the sheet to convert",
        )
        yield ConvArg(
            names=["--out-sep", ],
            default=cls.DEFAULT_OUT_SEP,
            help="The separator used in the output file",
        )
        yield ConvArg(
            names=["--line-terminator", ],
            default=cls.DEFAULT_LINE_TERMINATOR,
            help="The line terminator used in the output file",
        )
