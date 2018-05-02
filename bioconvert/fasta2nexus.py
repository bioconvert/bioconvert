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
"""FASTA2NEXUS conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['FASTA2NEXUS']


class FASTA2NEXUS(ConvBase):
    """Converts a sequence alignment in :term:`FASTA` format to :term:`NEXUS` format

    Conversion is based on Bio Python modules

    """
    _default_method = 'goalign'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("conda")
    def _method_goalign(self, threads=None, *args, **kwargs):
        """
        Convert fasta file in Nexus format using goalign tool.
        https://github.com/fredericlemoine/goalign

        The fasta file must be an alignemnt file, yhis mean all the sequences must
        have the same length (with the gap) otherwise an error will be raised
        """
        self.install_tool('goalign')
        cmd = 'goalign reformat nexus -i {infile} -o {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
