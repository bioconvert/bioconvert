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
"""Convert :term:`GENBANK` to :term:`EMBL` format"""

from bioconvert import ConvBase
from bioconvert.io.genbank import Genbank
from bioconvert.core.decorators import requires, requires_nothing
from bioconvert.core.decorators import  compressor

__all__ = ["GENBANK2FASTA"]


class GENBANK2FASTA(ConvBase):
    """Convert :term:`GENBANK` file to :term:`FASTA` file

    Methods are based on biopython [BIOPYTHON]_, squizz [SQUIZZ] and our
    own Bioconvert implementation.

    """
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GENBANK file
        :param str outfile: output EMBL filename

        """
        super(GENBANK2FASTA, self).__init__(infile, outfile, *args, **kargs)

        # squizz works as well but keeps lower cases while biopython uses upper
        # cases

    @requires("squizz")
    @compressor
    def _method_squizz(self, *args, **kwargs):
        """Header is less informative than the one obtained with biopython"""
        cmd = "squizz -f genbank -c fasta {} > {} ".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        from Bio import SeqIO
        SeqIO.convert(self.infile, "genbank", self.outfile, "fasta")

    @requires_nothing
    @compressor
    def _method_python(self, *args, **kwargs):
        reader = Genbank(self.infile)

        with open(self.outfile, "w") as writer:
            for idx, entry in enumerate(reader.read()):
                if "ORIGIN" in entry:
                    writer.write(">{} {}\n{}\n".format(
                        entry["VERSION"]["id"] if "VERSION" in entry else entry["LOCUS"]["id"],
                        entry["DEFINITION"] if "DEFINITION" in entry else "",
                        entry["ORIGIN"]))
                else:
                    print("Impossible to create a sequence for the entry number {}. Sequence not found after the keyword ORIGIN".format(idx))
