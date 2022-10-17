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
"""Convert :term:`EMBL` file to :term:`GENBANK` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import compressor, requires

__all__ = ["EMBL2GENBANK"]


class EMBL2GENBANK(ConvBase):
    """Convert :term:`EMBL` file to :term:`GENBANK` file


    Methods available are based on squizz [SQUIZZ]_ and
    biopython [BIOPYTHON]_.
    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input EMBL file
        :param str outfile: output GENBANK filename

        """
        super(EMBL2GENBANK, self).__init__(infile, outfile, *args, **kargs)

    @requires(external_binary="squizz")
    @compressor
    def _method_squizz(self, *args, **kwargs):
        """Header is less informative than the one obtained with biopython."""
        cmd = "squizz -f embl -c genbank {} > {} ".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """For this method we use the biopython package Bio.SeqIO.

        `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
        from Bio import SeqIO

        SeqIO.convert(self.infile, "embl", self.outfile, "genbank")
