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
"""convert :term:`CSV` to :term:`XLS` format"""
import csv

import colorlog

from bioconvert import ConvBase
from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import compressor, in_gz, requires, requires_nothing

logger = colorlog.getLogger(__name__)


__all__ = ["CSV2XLS"]


class CSV2XLS(ConvBase):
    """Convert :term:`CSV` file to :term:`XLS` file

    Methods available are based on python, pyexcel [PYEXCEL]_,
    or pandas [PANDAS]_.

    """

    #: Default value
    _default_method = "pandas"
    DEFAULT_IN_SEP = ","
    DEFAULT_LINE_TERMINATOR = "\n"
    DEFAULT_SHEET_NAME = "Sheet 1"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input CSV file
        :param str outfile: output XLS filename

        """
        super(CSV2XLS, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_libraries=["pyexcel", "pyexcel-xls"])
    @compressor
    def _method_pyexcel(self, in_sep=DEFAULT_IN_SEP, sheet_name=DEFAULT_SHEET_NAME, *args, **kwargs):
        """Do the conversion :term:`CSV` -> :term:`XLS` using pyexcel modules.

        `pyexcel documentation <http://docs.pyexcel.org/en/latest/>`_"""
        rows = []
        with open(self.infile, "r") as in_stream:
            reader = csv.reader(in_stream, delimiter=in_sep)
            for row in reader:
                rows.append(row)

        from collections import OrderedDict

        from pyexcel_xls import save_data

        data = OrderedDict()
        data.update({sheet_name: rows})
        save_data(self.outfile, data)

    @requires(python_libraries=["pandas"])
    @compressor
    def _method_pandas(self, in_sep=DEFAULT_IN_SEP, sheet_name=DEFAULT_SHEET_NAME, *args, **kwargs):
        """Do the conversion :term:`CSV` -> :term:`XLS` using Panda modules.

        `pandas documentation <https://pandas.pydata.org/docs/>`_"""
        import pandas as pd

        writer = pd.ExcelWriter(self.outfile, engine="openpyxl")
        pd.read_csv(self.infile, sep=in_sep, header="infer",).to_excel(
            excel_writer=writer,
            sheet_name=sheet_name,
            index=False,
        )
        writer.save()

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names=[
                "--sheet-name",
            ],
            default=cls.DEFAULT_SHEET_NAME,
            help="The name of the sheet to create",
        )
        yield ConvArg(
            names=[
                "--in-sep",
            ],
            default=cls.DEFAULT_IN_SEP,
            help="The separator used in the input file",
        )
