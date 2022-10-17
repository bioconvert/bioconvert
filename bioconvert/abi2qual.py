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
"""Convert :term:`ABI` format to :term:`QUAL` format"""
from bioconvert import ConvBase, requires

__all__ = ["ABI2QUAL"]


class ABI2QUAL(ConvBase):
    """Convert :term:`ABI` file to :term:`QUAL` file

    :term:`ABI` files are created by ABI sequencing machine and
    includes PHRED quality scores for base calls. This allows
    the creation of :term:`QUAL` files.

    Method implemented is based on BioPython [BIOPYTHON]_.

    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input ABI file
        :param str outfile: output QUAL filename

        """
        super(ABI2QUAL, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        """For this method we use the biopython package Bio.SeqIO.

        `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
        from Bio import SeqIO

        records = SeqIO.parse(self.infile, "abi")
        # output using SeqIO.write(records, self.outfile, "qual") is not
        # standard so we write our own conversion here below
        with open(self.outfile, "w") as fout:
            for rec in records:
                header = rec.name
                qual = rec.letter_annotations["phred_quality"]
                qual = "".join([str(x) for x in qual])
                fout.write(">{}\n".format(header))
                fout.write("{}\n".format(qual))
