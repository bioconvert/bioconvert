# -*- coding: utf-8 -*-
##############################################################################
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
"""convert :term:`CSV` to :term:`XLS` format"""
import csv

import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires, requires_nothing
from bioconvert.core.decorators import compressor, in_gz
from bioconvert.core.base import ConvArg

logger = colorlog.getLogger(__name__)


__all__ = ["CSV2XLS"]


class CSV2XLS(ConvBase):
    """Convert :term:`CSV` file to :term:`XLS` file

    Methods available are based on python, pyexcel [PYEXCEL]_,
    or pandas [PANDAS]_.

    """
    _default_method = "pandas"
    DEFAULT_IN_SEP = ','
    DEFAULT_LINE_TERMINATOR = '\n'
    DEFAULT_SHEET_NAME = "Sheet 1"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input CSV file
        :param str outfile: output XLS filename

        """
        super(CSV2XLS, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_libraries=["pyexcel", "pyexcel-xls"])
    @compressor
    def _method_pyexcel(
            self,
            in_sep=DEFAULT_IN_SEP,
            sheet_name=DEFAULT_SHEET_NAME,
            *args, **kwargs):
        """
        Do the conversion :term:`CSV` -> :term:`XLS` using pyexcel modules

        """
        rows = []
        with open(self.infile, "r") as in_stream:
            reader = csv.reader(in_stream, delimiter=in_sep)
            for row in reader:
                rows.append(row)

        from pyexcel_xls import save_data
        from collections import OrderedDict
        data = OrderedDict()
        data.update({sheet_name: rows})
        save_data(self.outfile, data)

    @requires(python_libraries=["pandas"])
    @compressor
    def _method_pandas(
            self,
            in_sep=DEFAULT_IN_SEP,
            sheet_name=DEFAULT_SHEET_NAME,
            *args, **kwargs):
        """
        Do the conversion :term:`CSV` -> :term:`XLS` using Panda modules

        """
        import pandas as pd
        writer = pd.ExcelWriter(self.outfile)
        pd.read_csv(
            self.infile,
            sep=in_sep,
            header='infer',
        ) \
            .to_excel(
            excel_writer=writer,
            sheet_name=sheet_name,
            index=False,
        )
        writer.save()

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names=["--sheet-name", ],
            default=cls.DEFAULT_SHEET_NAME,
            help="The name of the sheet to create",
        )
        yield ConvArg(
            names=["--in-sep", ],
            default=cls.DEFAULT_IN_SEP,
            help="The separator used in the input file",
        )
