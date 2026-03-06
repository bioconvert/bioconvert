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
"""Convert :term:`GENBANK` to :term:`FAA` format"""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import compressor, requires

_log = colorlog.getLogger(__name__)

__all__ = ["GENBANK2FAA"]


class GENBANK2FAA(ConvBase):
    """Convert :term:`GENBANK` file to :term:`FAA` file

    Extract protein sequences from CDS features with a ``/translation``
    qualifier and write them in :term:`FASTA` amino acid format.

    Methods available are based on biopython [BIOPYTHON]_.

    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GENBANK file
        :param str outfile: output FAA filename
        """
        super(GENBANK2FAA, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """Extract CDS translations using :class:`Bio.SeqIO`.

        For each CDS feature that contains a ``/translation`` qualifier the
        protein identifier and product name are used as the FASTA header and
        the translated sequence is written to the output file.

        `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
        from Bio import SeqIO

        with open(self.outfile, "w") as fout:
            for record in SeqIO.parse(self.infile, "genbank"):
                for feature in record.features:
                    if feature.type == "CDS":
                        qualifiers = feature.qualifiers
                        if "translation" not in qualifiers:
                            continue
                        protein_id = qualifiers.get("protein_id", ["unknown"])[0]
                        product = qualifiers.get("product", [""])[0]
                        translation = qualifiers["translation"][0]
                        header = protein_id
                        if product:
                            header = "{} {}".format(protein_id, product)
                        fout.write(">{}\n{}\n".format(header, translation))
