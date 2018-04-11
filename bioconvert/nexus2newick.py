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
"""NEXUS2NEWICK conversion"""
import os

from bioconvert import ConvBase

import colorlog

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['NEXUS2NEWICK']


class NEXUS2NEWICK(ConvBase):
    """
    Converts a tree file from :term:`NEXUS` format to :term:`NEWICK` format. ::
    """
    _default_method = 'gotree'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEXUS` file.
        :param str outfile: (optional) output :term:`NEWICK` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        _log.warning("biopython methods rounds up values (5 digits)")
        from Bio import Phylo
        Phylo.convert(self.infile, "nexus", self.outfile, "newick")

    @requires("conda")
    def _method_gotree(self, threads=None, *args, **kwargs):
        """
        Convert :term:`NEXUS`  file in :term:`NEWICK` format using gotree tool.
        https://github.com/fredericlemoine/gotree

        :param threads: not used.
        """
        self.install_tool('gotree')
        cmd = 'gotree reformat newick -i {infile} -o {outfile} -f nexus'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
