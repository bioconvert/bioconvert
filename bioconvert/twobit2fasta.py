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
"""Conversion from :term:`TWOBIT` to :term:`FASTA` format"""
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ['TWOBIT2FASTA']


class TWOBIT2FASTA(ConvBase):
    """Converts a sequence alignment in :term:`TWOBIT` format to :term:`FASTA` format

    Conversion is based on UCSC [UCSC]_ and py2bit.

    """
    _default_method = 'py2bit'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`TWOBIT` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        super(TWOBIT2FASTA, self).__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("twoBitToFa")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert twobit file in fasta format using ucsc twobittofa.
        https://genome.ucsc.edu/goldenPath/help/twoBit.html
        """
        cmd = 'twoBitToFa {infile} {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)

    # py2bit is from deeptols repo
    @requires(python_library="py2bit")
    def _method_py2bit(self, *args, **kwargs):

        import py2bit
        data = py2bit.open(self.infile)
        with open(self.outfile, "w") as fout:
            for chrom in sorted(data.chroms()):
                seq = data.sequence(chrom)
                fout.write(">{}\n{}\n".format(chrom, seq))

