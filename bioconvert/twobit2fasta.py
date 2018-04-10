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
"""TWOBIT2FASTA conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['TWOBIT2FASTA']


class TWOBIT2FASTA(ConvBase):
    """Converts a sequence alignment in :term:`TWOBIT` format to :term:`FASTA` format

    Conversion is based on UCSC twobit2fa

    """
    _default_method = 'ucsc'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`TWOBIT` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("twoBitToFa")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert twobit file in fasta format using ucsc twobittofa.
        https://genome.ucsc.edu/goldenPath/help/twoBit.html
        """
        cmd = 'twoBitToFa {infile} {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
