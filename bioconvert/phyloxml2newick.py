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
"""Converts :term:`PHYLOXML` file to :term:`NEWICK` format."""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.decorators import compressor

_log = colorlog.getLogger(__name__)


__all__ = ['PHYLOXML2NEWICK']


class PHYLOXML2NEWICK(ConvBase):
    """
    Converts a tree file from :term:`PHYLOXML` format to :term:`NEWICK` format. 

    Methods available are based on gotree [GOTREE]_.

    """
    _default_method = 'gotree'

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLOXML` file.
        :param str outfile: (optional) output :term:`NEWICK` file
        """
        super(PHYLOXML2NEWICK, self).__init__(infile, outfile)

    @requires("go")
    @compressor
    def _method_gotree(self, *args, **kwargs):
        """
        Convert :term:`PHYLOXML`  file in :term:`NEWICK` format using gotree tool.
        https://github.com/fredericlemoine/gotree

        """
        self.install_tool('gotree')
        cmd = 'gotree reformat newick -i {infile} -o {outfile} -f phyloxml'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
