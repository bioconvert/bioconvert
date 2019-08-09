# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
"""Convert :term:`BEDGRAPH` format to :term:`WIGGLE` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires
import colorlog

_log = colorlog.getLogger(__name__)

__all__ = ["BEDGRAPH2WIGGLE"]


class BEDGRAPH2WIGGLE(ConvBase):
    """Convert sorted :term:`BEDGRAPH` file into :term:`WIGGLE` file

    Methods available are based on wiggletools [WIGGLETOOLS]_.

    """
    _default_method = "wiggletools"

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BEDGRAPH file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super(BEDGRAPH2WIGGLE, self).__init__(infile, outfile)

    @requires("wiggletools")
    def _method_wiggletools(self, *args, **kwargs):
        """wiggletools based method

        Extension must be .bg
        """
        if self.infile.endswith(".bg") is False:
            raise ValueError("The input bedgraph file must have the .bg extension (wiggletools requirements)")
        cmd = "wiggletools {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

