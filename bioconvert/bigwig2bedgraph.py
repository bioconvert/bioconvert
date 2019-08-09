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
"""Convert :term:`BIGWIG` to :term:`BEDGRAPH` format """
import os

import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BIGWIG2BEDGRAPH"]


class BIGWIG2BEDGRAPH(ConvBase):
    """Converts a sequence alignment in :term:`BIGWIG` format to :term:`BEDGRAPH` format

    Conversion is based on ucsc bigWigToBedGraph tool or pybigwig (default)
    [DEEPTOOLS]_.

    """
    _default_method = 'pybigwig'

    def __init__(self, infile, outfile):#=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`BIGWIG` file.
        :param str outfile: (optional) output :term:`BEDGRAPH` file
        """
        super(BIGWIG2BEDGRAPH, self).__init__(infile, outfile)
        #self.alphabet = alphabet

    @requires("bigWigToBedGraph")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert bigwig file in bedgraph format using ucsc tool.
        https://genome.ucsc.edu/goldenPath/help/bedgraph.html
        """
        cmd = 'bigWigToBedGraph {infile}  {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)

    @requires(python_library="pyBigWig")
    def _method_pybigwig(self, *args, **kwargs):
        import pyBigWig
        bw = pyBigWig.open(self.infile)
        assert bw.isBigWig() is True, "Not a valid bigWig file"

        with open(self.outfile, "w") as fout:
            for chrom in  bw.chroms():
                for tup in bw.intervals(chrom):
                    s, e, val = tup
                    if int(val) == val:
                        val = int(val)
                    fout.write("{}\t{}\t{}\t{}\n".format(chrom, s, e, val))



