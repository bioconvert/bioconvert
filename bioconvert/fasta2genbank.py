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
"""Convert :term:`FASTA` to :term:`GENBANK` format"""
import datetime
from math import log, floor

from bioconvert import ConvBase
from bioconvert.core.decorators import requires, requires_nothing
from bioconvert.core.decorators import compressor
from bioconvert.io.fasta import Fasta


__all__ = ["FASTA2GENBANK"]


class FASTA2GENBANK(ConvBase):
    """Convert :term:`FASTA` file to :term:`GENBANK` file

    Methods available are based on squizz [SQUIZZ]_ or biopython [BIOPYTHON]_ or 
    Bioconvert pure implementation (default).

    """

    # squizz works as well but keeps lower cases while 
    # biopython uses upper cases
    _default_method = "bioconvert"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input FASTA file
        :param str outfile: output GENBANK filename

        """
        super(FASTA2GENBANK, self).__init__(infile, outfile, *args, **kargs)

    @requires("squizz")
    @compressor
    def _method_squizz(self, *args, **kwargs):
        """Header is less informative than the one obtained with biopython"""
        cmd = "squizz -f fasta -c genbank  {} > {} ".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        print("Using DNA alphabet for now")
        from Bio import SeqIO, Alphabet
        SeqIO.convert(self.infile, "fasta", self.outfile, "genbank",
            alphabet=Alphabet.generic_dna)

    # --- Pure python methods ---

    @requires_nothing
    def _method_bioconvert(self, *args, **kwargs):
        reader = Fasta(self.infile)

        with open(self.outfile, "w") as writer:
            for sequence in reader.read():
                seq_size = len(sequence["value"])
                num_digit = floor(log(seq_size, 10)) + 1

                # Sequence header
                now = datetime.datetime.now()
                writer.write("LOCUS       {}{}{} bp XXXXXX              XXX {}-{}-{}\n".format(
                    sequence["id"],
                    " "*(max(1, 28 - len(sequence["id"]) - num_digit)),
                    seq_size,
                    now.day, now.month, now.year))
                writer.write("DEFINITION  {}\n".format(sequence["comment"]))
                writer.write("ORIGIN      \n")

                # Print sequence
                for seq_idx in range(0, seq_size, 60):
                    # Write line header (idx in the sequence)
                    idx_num_digit = floor(log(seq_idx+1, 10)) + 1
                    writer.write("{}{}".format(" "*(9 - idx_num_digit), seq_idx+1))

                    # write the sequence itself
                    for i in range(6):
                        begin = seq_idx+i*10
                        end = seq_idx+(i+1)*10

                        # sequence over before this slice
                        if begin >= seq_size:
                            break
                        # sequence over during this slice
                        elif end > seq_size:
                            writer.write(" {}".format(sequence["value"][begin:seq_size]))
                        else:
                            writer.write(" {}".format(sequence["value"][begin:end]))

                    # newline
                    writer.write("\n")
                writer.write("//\n")
