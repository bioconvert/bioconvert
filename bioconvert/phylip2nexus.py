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
"""PHYLIP2NEXUS conversion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['PHYLIP2NEXUS']


class PHYLIP2NEXUS(ConvBase):
    """
    Converts a sequence alignment from :term:`PHYLIP` format to :term:`NEXUS` format. ::
    """
    _default_method = 'goalign'

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super().__init__(infile, outfile)

    @requires("conda")
    def _method_goalign(self, *args, **kwargs):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`NEXUS` format using goalign tool.
        https://github.com/fredericlemoine/goalign

        """
        self.install_tool('goalign')
        cmd = 'goalign reformat nexus -i {infile} -o {outfile} -p'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
