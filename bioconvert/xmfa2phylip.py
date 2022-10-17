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
"""Convert :term:`XMFA` to :term:`PHYLIP` format"""

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase, compressor, requires

_log = colorlog.getLogger(__name__)


__all__ = ["XMFA2PHYLIP"]


class XMFA2PHYLIP(ConvBase):
    """Converts a sequence alignment from :term:`XMFA` to :term:`PHYLIP` format.

    Method available based on biopython [BIOPYTHON]_.

    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`PHYLIP` file.
        :param str outfile: (optional) output :term:`NEXUS` file
        """
        super(XMFA2PHYLIP, self).__init__(infile, outfile)

    @requires(python_libraries=["biopython"])
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """Convert :term:`XMFA` interleaved file in :term:`PHYLIP` (Mauve)format.

        `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
        sequences = list(SeqIO.parse(self.infile, "mauve"))
        count = SeqIO.write(sequences, self.outfile, "phylip")
        _log.info("Converted %d records to xmfa" % count)
