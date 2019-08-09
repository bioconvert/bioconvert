# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
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
###########################################################################
"""Convert :term:`XLS` format to :term:`CSV` format"""
import csv

import colorlog

from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires, requires_nothing

from bioconvert import ConvBase

logger = colorlog.getLogger(__name__)


class ODS2CSV(ConvBase):
    """Convert :term:`XLS` file into :term:`CSV` file

    Method based on pyexcel [PYEXCEL].

    """
    _default_method = "pyexcel"
    DEFAULT_OUT_SEP = ','
    DEFAULT_LINE_TERMINATOR = '\n'

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:
        """
        super(ODS2CSV, self).__init__(infile, outfile)

    @requires(python_libraries=["pyexcel", "pyexcel-ods3"])
    def _method_pyexcel(
            self,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            sheet_name=0,
            *args, **kwargs):
        """Do the conversion :term:`XLS` -> :term:`CSV` using Panda library"""

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
