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
"""Convert :term:`CSV` format to :term:`TSV` file"""

import colorlog

from bioconvert import tsv2csv
from bioconvert.core.decorators import requires, requires_nothing
from bioconvert.core.decorators import compressor, in_gz


logger = colorlog.getLogger(__name__)


class CSV2TSV(tsv2csv.TSV2CSV):
    """Convert :term:`CSV` file into :term:`TSV` file

    Available methods: Python, Pandas

    .. plot::

        from bioconvert.csv2tsv import CSV2TSV
        from bioconvert import bioconvert_data, logger
        from easydev import TempFile

        logger.level = 'CRITICAL'
        with TempFile(suffix=".csv") as fh:
           infile = bioconvert_data("test_tabulated.tsv")
           convert = CSV2TSV(infile, fh.name)
           convert.boxplot_benchmark(N=50)

    .. seealso:: :class:`~bioconvert.csv2tsv.TSV2CSV`
    """
    _default_method = "python"
    DEFAULT_IN_SEP = ','
    DEFAULT_OUT_SEP = '\t'
    DEFAULT_LINE_TERMINATOR = '\n'

    def __init__(self, infile, outfile):
        """.. rubric:: Constructor

        :param str infile:
        :param str outfile:
        """
        super().__init__(infile, outfile)

    @requires_nothing
    @compressor
    def _method_python(
            self,
            in_sep=DEFAULT_IN_SEP,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            *args, **kwargs):
        """
        Do the conversion :term:`CSV` -> :term:`TSV` using standard Python modules
        """
        super()._method_python(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)

    @requires_nothing
    @compressor
    def _method_python_v2(
            self,
            in_sep=DEFAULT_IN_SEP,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            *args, **kwargs):
        """
        Do the conversion :term:`CSV` -> :term:`CSV` using csv module

        .. note:: This method cannot escape nor quote output char
        """
        super()._method_python_v2(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)

    @requires(python_library="pandas")
    @compressor
    def _method_pandas(
            self,
            in_sep=DEFAULT_IN_SEP,
            out_sep=DEFAULT_OUT_SEP,
            line_terminator=DEFAULT_LINE_TERMINATOR,
            *args, **kwargs):
        """
        Do the conversion :term:`CSV` -> :term:`TSV` using Pandas library
        """
        super()._method_pandas(in_sep=in_sep, out_sep=out_sep, *args, **kwargs)
