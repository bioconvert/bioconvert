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
"""NEXUS2PHYLOXML conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, extensions

_log = colorlog.getLogger(__name__)


__all__ = ['NEXUS2PHYLOXML']


class NEXUS2PHYLOXML(ConvBase):
    """
    Converts a tree file from :term:`NEXUS` format to :term:`PHYLOXML` format. ::
    """


    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEXUS` file.
        :param str outfile: (optional) output :term:`PHYLOXML` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet
        self._default_method = 'gotree'

    def _method_gotree(self, threads=None, *args, **kwargs):
        """
        Convert :term:`NEXUS`  file in :term:`PHYLOXML` format using gotree tool.
        https://github.com/fredericlemoine/gotree

        :param threads: not used.
        """
        self.install_tool('gotree')
        cmd = 'gotree reformat phyloxml -i {infile} -o {outfile} -f nexus'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
