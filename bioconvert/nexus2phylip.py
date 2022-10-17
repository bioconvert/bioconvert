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
"""Converts :term:`NEXUS` file to :term:`PHYLIP` format."""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import compressor, requires

_log = colorlog.getLogger(__name__)


__all__ = ["NEXUS2PHYLIP"]


class NEXUS2PHYLIP(ConvBase):
    """
    Converts a sequence alignment from :term:`NEXUS` format to :term:`PHYLIP` format.

    Methods available are based on goalign [GOALIGN]_.

    """

    #: Default value
    _default_method = "goalign"

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEXUS` file.
        :param str outfile: (optional) output :term:`PHYLIP` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("go")
    @compressor
    def _method_goalign(self, *args, **kwargs):
        """Convert :term:`NEXUS` interleaved file in :term:`PHYLIP` format using goalign tool.

        `goalign documentation <https://github.com/fredericlemoine/goalign>`_"""
        self.install_tool("goalign")
        cmd = "goalign reformat phylip -i {infile} -o {outfile} -x".format(infile=self.infile, outfile=self.outfile)
        self.execute(cmd)
