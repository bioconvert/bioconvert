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
"""Convert :term:`CLUSTAL` to :term:`NEXUS` format"""
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import compressor, requires

_log = colorlog.getLogger(__name__)

__all__ = ["CLUSTAL2NEXUS"]


class CLUSTAL2NEXUS(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`NEXUS` format. ::

    Methods available are based on squizz [SQUIZZ] or biopython [BIOPYTHON], and
    goalign [GOALIGN].

    """

    #: Default value
    _default_method = "goalign"

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super(CLUSTAL2NEXUS, self).__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("go")
    @compressor
    def _method_goalign(self, threads=None, *args, **kwargs):
        """Convert :term:`CLUSTAL` file in  :term:`NEXUS` format using goalign tool.

        `goalign docuMentation <https://github.com/fredericlemoine/goalign>`_"""
        self.install_tool("goalign")
        cmd = "goalign reformat nexus --clustal -i {infile} -o {outfile}".format(
            infile=self.infile, outfile=self.outfile
        )
        self.execute(cmd)
