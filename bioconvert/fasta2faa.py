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
# You should have received a copy of the GNT General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################
"""Convert :term:`FASTA` format to :term:`FAA` format"""
import textwrap
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import compressor
from bioconvert.core.decorators import requires
from bioconvert.io.fasta import Fasta

_log = colorlog.getLogger(__name__)


__all__ = ["FASTA2FAA"]


class FASTA2FAA(ConvBase):
    """

    Methods available is a bioconvert implementation.

    """
    _default_method = "bioconvert"

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FASTA file
        :param str outfile: The path to the output FASTQ file
        """
        super(FASTA2FAA, self).__init__(infile, outfile)

        self.codons = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
            'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        }


    @compressor
    def _method_bioconvert(self, *args, **kwargs):
        data = Fasta(self.infile)
        with open(self.outfile, "w") as fout:
            for line in data.read():
                fout.write(">{}\t{}\n".format(line['id'], line['comment']))
                seq = line['value']

                # translate the sequence
                aa = "".join((self.codons.get(seq[i:i+3], "X") for i in range(0, len(seq), 3)))
                fout.write("\n".join(textwrap.wrap(aa, 60)) + "\n")




