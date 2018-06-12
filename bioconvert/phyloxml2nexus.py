# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################

"""PHYLOXML2NEXUS convertion"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['PHYLOXML2NEXUS']


class PHYLOXML2NEXUS(ConvBase):
    """
    Converts a tree file from :term:`PHYLOXML` format to :term:`NEXUS` format. ::
    """
    _default_method = 'gotree'


    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLOXML` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("conda")
    def _method_gotree(self, threads=None, *args, **kwargs):
        """
        Convert :term:`PHYLOXML`  file in :term:`NEXUS` format using gotree tool.
        https://github.com/fredericlemoine/gotree

        :param threads: not used.
        """
        self.install_tool('gotree')
        cmd = 'gotree reformat nexus -i {infile} -o {outfile} -f phyloxml'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
