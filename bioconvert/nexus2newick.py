###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""Converts :term:`NEXUS` file to :term:`NEWICK` format."""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import compressor, requires

_log = colorlog.getLogger(__name__)


__all__ = ["NEXUS2NEWICK"]


class NEXUS2NEWICK(ConvBase):
    """
    Converts a tree file from :term:`NEXUS` format to :term:`NEWICK` format.

    Methods available are based on biopython [BIOPYTHON]_ or
    goalign [GOALIGN]_.

    """

    #: Default value
    _default_method = "gotree"

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEXUS` file.
        :param str outfile: (optional) output :term:`NEWICK` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """For this method we use the biopython package Bio.Phylo.

        `Bio.Phylo Documentation <https://biopython.org/wiki/Phylo>`_"""
        _log.warning("biopython methods rounds up values (5 digits)")
        from Bio import Phylo

        Phylo.convert(self.infile, "nexus", self.outfile, "newick")

    @requires("gotree")
    @compressor
    def _method_gotree(self, *args, **kwargs):
        """Convert :term:`NEXUS`  file in :term:`NEWICK` format using gotree tool.

        `gotree documentation <https://github.com/fredericlemoine/gotree>`_"""
        self.install_tool("gotree")
        cmd = "gotree reformat newick -i {infile} -o {outfile} -f nexus".format(
            infile=self.infile, outfile=self.outfile
        )
        self.execute(cmd)
