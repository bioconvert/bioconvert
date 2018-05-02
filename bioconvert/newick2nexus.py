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
"""NEWICKTONEXUS conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['NEWICK2NEXUS']


class NEWICK2NEXUS(ConvBase):
    """
    Converts a tree file from :term:`NEWICK` format to :term:`NEXUS` format. ::
    """
    _default_method = 'gotree'

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEWICK` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super().__init__(infile, outfile)

    @requires("conda")
    def _method_gotree(self, threads=None, *args, **kwargs):
        """
        Convert :term:`NEWICK`  file in :term:`NEXUS` format using gotree tool.
        https://github.com/fredericlemoine/gotree

        :param threads: not used.
        """
        self.install_tool('gotree')
        cmd = 'gotree reformat nexus -i {infile} -o {outfile} -f newick'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
