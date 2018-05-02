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
""" description """
import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)

__all__ = ["CLUSTAL2PHYLIP"]


class CLUSTAL2PHYLIP(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`PHYLIP` format. ::

        converter = CLUSTAL2PHYLIP(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`PHYLIP` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    def _method_biopython(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using biopython.

        :param threads: not used.
        """
        sequences = list(SeqIO.parse(self.infile, "clustal", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "phylip")
        _log.info("Converted %d records to phylip" % count)

    @requires("squizz")
    def _method_squizz(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c PHYLIPI {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)

