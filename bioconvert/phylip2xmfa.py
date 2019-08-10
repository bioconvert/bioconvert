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
"""Converts :term:`PHYLIP` file to :term:`XMFA` format."""
import colorlog
from Bio import AlignIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.decorators import compressor

_log = colorlog.getLogger(__name__)


__all__ = ['PHYLIP2XMFA']


class PHYLIP2XMFA(ConvBase):
    """
    Converts a sequence alignment from :term:`PHYLIP` format to :term:`XMFA` 

    Methods available are based on biopython [BIOPYTHON]_.

    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super().__init__(infile, outfile)

    @requires(python_libraries=["biopython"])
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """
        Convert :term:`PHYLIP` interleaved file in :term:`XMFA` (Mauve)format.

        """
        sequences = list(AlignIO.parse(self.infile, "phylip"))
        count = AlignIO.write(sequences, self.outfile, "mauve")
        _log.info("Converted %d records to phylip" % count)


