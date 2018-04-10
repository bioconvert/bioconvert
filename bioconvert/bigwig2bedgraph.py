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
"""BIGWIG2BEDGRAPH conversion """
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BIGWIG2BEDGRAPH"]


class BIGWIG2BEDGRAPH(ConvBase):
    """Converts a sequence alignment in :term:`BIGWIG` format to :term:`BEDGRAPH` format

    Conversion is based on ucsc bigWigToBedGraph tool

    """
    _default_method = 'ucsc'

    def __init__(self, infile, outfile):#=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`BIGWIG` file.
        :param str outfile: (optional) output :term:`BEDGRAPH` file
        """
        super().__init__(infile, outfile)
        #self.alphabet = alphabet

    @requires("bigWigToBedGraph")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert bigwig file in bedgraph format using ucsc tool.
        https://genome.ucsc.edu/goldenPath/help/bedgraph.html
        """
        cmd = 'bigWigToBedGraph {infile}  {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
