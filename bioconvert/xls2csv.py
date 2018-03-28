"""Convert :term:`XLS` format to :term:`CSV` file"""
import csv

import colorlog

from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires, requires_nothing

from bioconvert import ConvBase

logger = colorlog.getLogger(__name__)


class XLS2CSV(ConvBase):
    """Convert :term:`XLS` file into :term:`CSV` file

    """
    _default_method = "pandas"
    DEFAULT_OUT_SEP = ','
    DEFAULT_LINE_TERMINATOR = '\n'

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        """
        super().__init__(infile, outfile)

    @requires(python_libraries=["pandas", "xlrd"])
    def _method_panda(
            self,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            sheet_name=0,
            *args, **kwargs):
        """
        do the conversion :term`XLS` -> :term:'CSV` using Panda modules

        """
        import pandas as pd
        pd.read_excel(
            self.infile,
            sheet_name=sheet_name,
        ) \
            .to_csv(
            self.outfile,
            sep=out_sep,
            line_terminator=line_terminator,
            index=False,
            header='infer'
        )

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
