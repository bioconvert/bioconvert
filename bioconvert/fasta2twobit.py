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
"""FASTA2TWOBIT conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['FASTA2TWOBIT']


class FASTA2TWOBIT(ConvBase):
    """Converts a sequence alignment in :term:`FASTA` format to :term:`TWOBIT` format

    Conversion is based on UCSC faToTwoBit

    """
    _default_method = 'ucsc'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`TWOBIT` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("faToTwoBit")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert fasta file in twobit format using ucsc faToTwoBit.
        https://genome.ucsc.edu/goldenPath/help/twoBit.html
        """
        cmd = 'faToTwoBit {infile} {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
