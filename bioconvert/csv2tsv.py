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
"""Convert :term:`CSV` format to :term:`TSV` format"""
import csv

import colorlog

from bioconvert.core.base import ConvArg, ConvBase
from bioconvert.core.decorators import compressor, in_gz, requires, requires_nothing

logger = colorlog.getLogger(__name__)


class CSV2TSV(ConvBase):
    """Convert :term:`CSV` file into :term:`TSV` file

    Available methods: Python, Pandas

    Methods available are based on python or Pandas [PANDAS]_.

    .. seealso:: :class:`~bioconvert.csv2tsv.TSV2CSV`
    """

    #: Default value
    _default_method = "python"
    DEFAULT_IN_SEP = ","
    DEFAULT_OUT_SEP = "\t"
    DEFAULT_LINE_TERMINATOR = "\n"

    def __init__(self, infile, outfile):
        """.. rubric:: Constructor

        :param str infile: comma-separated file
        :param str outfile: tabulated file
        """
        super(CSV2TSV, self).__init__(infile, outfile)

    @requires_nothing
    @compressor
    def _method_python(
        self, in_sep=DEFAULT_IN_SEP, out_sep=DEFAULT_OUT_SEP, line_terminator=DEFAULT_LINE_TERMINATOR, *args, **kwargs
    ):
        """Do the conversion :term:`CSV` -> :term:`TSV` using standard Python modules.

        `csv documentation <https://docs.python.org/3/library/csv.html>`_"""
        with open(self.infile, "r") as in_stream, open(self.outfile, "w") as out_stream:
            writer = csv.writer(out_stream, delimiter=out_sep, lineterminator=line_terminator)
            reader = csv.reader(in_stream, delimiter=in_sep)
            for row in reader:
                writer.writerow(row)

    @requires_nothing
    @compressor
    def _method_python_v2(
        self, in_sep=DEFAULT_IN_SEP, out_sep=DEFAULT_OUT_SEP, line_terminator=DEFAULT_LINE_TERMINATOR, *args, **kwargs
    ):
        """Do the conversion :term:`CSV` -> :term:`CSV` using csv module.

        .. note:: This method cannot escape nor quote output char

        `csv documentation <https://docs.python.org/3/library/csv.html>`_"""
        with open(self.infile, "r") as in_stream, open(self.outfile, "w") as out_stream:
            reader = csv.reader(in_stream, delimiter=in_sep)
            for row in reader:
                out_stream.write(out_sep.join(row))
                out_stream.write(line_terminator)

    @requires(python_library="pandas")
    @compressor
    def _method_pandas(
        self, in_sep=DEFAULT_IN_SEP, out_sep=DEFAULT_OUT_SEP, line_terminator=DEFAULT_LINE_TERMINATOR, *args, **kwargs
    ):
        """Do the conversion :term:`CSV` -> :term:`TSV` using Pandas library


        `pandas documentation <https://pandas.pydata.org/docs/>`_"""
        import pandas as pd

        pd.read_csv(self.infile, sep=in_sep,).to_csv(
            self.outfile,
            sep=out_sep,
            line_terminator=line_terminator,
            index=False,
            header="infer",
        )

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names=[
                "--in-sep",
            ],
            default=cls.DEFAULT_IN_SEP,
            help="The separator used in the input file",
        )
        yield ConvArg(
            names=[
                "--out-sep",
            ],
            default=cls.DEFAULT_OUT_SEP,
            help="The separator used in the output file",
        )
        yield ConvArg(
            names=[
                "--line-terminator",
            ],
            default=cls.DEFAULT_LINE_TERMINATOR,
            help="The line terminator used in the output file",
        )
