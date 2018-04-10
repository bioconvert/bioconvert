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
"""PHYLIP2FASTA conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['PHYLIP2FASTA']


class PHYLIP2FASTA(ConvBase):
    """Converts a sequence alignment in :term:`PHYLIP` format to :term:`FASTA` format

    Conversion is based on Bio Python modules

    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        sequences = list(SeqIO.parse(self.infile, "phylip", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "fasta")
        _log.debug("Converted %d records to fasta" % count)

    @requires("squizz")
    def _method_squizz(self, *args, **kwargs):
        """
        Convert Phylip inteleaved file in fasta format using squizz tool.
        The fasta file is an alignemnt, that means the gap are conserved.
        """
        cmd = 'squizz -c FASTA {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)

    @requires("conda")
    def _method_goalign(self, *args, **kwargs):
        """
        Convert fasta file in Phylip interleaved format using goalign tool.
        https://github.com/fredericlemoine/goalign
        The fasta file must be an alignemnt file, yhis mean all the sequences must
        have the same length (with the gap) otherwise an error will be raised
        """
        self.install_tool('goalign')
        cmd = 'goalign reformat fasta -i {infile} -p -o {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
