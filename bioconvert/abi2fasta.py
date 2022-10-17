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
"""Convert :term:`ABI` format to :term:`FASTA` format"""
from bioconvert import ConvBase, requires

__all__ = ["ABI2FASTA"]


class ABI2FASTA(ConvBase):
    """Convert :term:`ABI` file to :term:`FASTQ` file

    :term:`ABI` files are created by ABI sequencing machine and includes
    PHRED quality scores for base calls. This allows the creation of
    :term:`FastA` files.

    Method implemented is based on BioPython [BIOPYTHON]_.

    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input ABI file
        :param str outfile: output FASTA filename

        """
        super(ABI2FASTA, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        """For this method we use the biopython package Bio.SeqIO.

        :reference: `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
        from Bio import SeqIO

        records = SeqIO.parse(self.infile, "abi")
        SeqIO.write(records, self.outfile, "fasta")
