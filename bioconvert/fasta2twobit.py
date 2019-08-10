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
"""Convert :term:`FASTA` to :term:`TWOBIT` format"""
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.decorators import compressor

_log = colorlog.getLogger(__name__)


__all__ = ['FASTA2TWOBIT']


class FASTA2TWOBIT(ConvBase):
    """Converts a sequence alignment in :term:`FASTA` format to :term:`TWOBIT` format

    Methods available are based on UCSC faToTwoBit [UCSC]_.

    """
    _default_method = 'ucsc'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`TWOBIT` file
        """
        super(FASTA2TWOBIT, self).__init__(infile, outfile)

    @requires("faToTwoBit")
    @compressor
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert fasta file in twobit format using ucsc faToTwoBit.
        https://genome.ucsc.edu/goldenPath/help/twoBit.html
        """
        cmd = 'faToTwoBit {infile} {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)



