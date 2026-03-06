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
"""Convert :term:`TRANSFAC` to :term:`JASPAR` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

__all__ = ["TRANSFAC2JASPAR"]


class TRANSFAC2JASPAR(ConvBase):
    """Convert :term:`TRANSFAC` motif file to :term:`JASPAR` motif file

    Conversion is based on the biopython [BIOPYTHON]_ library.

    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`TRANSFAC` file
        :param str outfile: output :term:`JASPAR` file
        """
        super(TRANSFAC2JASPAR, self).__init__(infile, outfile, *args, **kwargs)

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        """Convert :term:`TRANSFAC` motif file to :term:`JASPAR` format using biopython.

        `Bio.motifs Documentation <https://biopython.org/docs/latest/api/Bio.motifs.html>`_"""
        from Bio import motifs

        with open(self.infile, "r") as infile:
            record = motifs.parse(infile, "transfac")

        # Preserve AC (accession) and ID (name) metadata from TRANSFAC fields
        for motif in record:
            if not motif.name:
                motif.name = motif.get("ID", "")
            if not getattr(motif, "matrix_id", None):
                motif.matrix_id = motif.get("AC", None)

        with open(self.outfile, "w") as outfile:
            outfile.write(motifs.write(record, "jaspar"))
