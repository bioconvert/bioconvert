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


__all__ = ['PHYLIP2STOCKHOLM']


class PHYLIP2STOCKHOLM(ConvBase):
    """
    Converts a sequence alignment from :term:`PHYLIP` interleaved format to :term:`STOCKHOLM` format. ::

        converter = PHYLIP2STOCKHOLM(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz
    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`STOCKHOLM` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    def _method_biopython(self, threads=None, *args, **kwargs):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`STOCKHOLM` format using biopython.

        :param threads: not used.
        """
        sequences = list(SeqIO.parse(self.infile, "phylip", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "stockholm")
        _log.info("Converted %d records to stockholm" % count)

    @requires("squizz")
    def _method_squizz(self, threads=None, *args, **kwargs):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`STOCKHOLM` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c STOCKHOLM {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
