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
"""Convert :term:`BEDGRAPH` to :term:`BIGWIG` format"""
import os
import colorlog

from bioconvert.core.base import ConvArg, ConvBase

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BEDGRAPH2BIGWIG"]


class BEDGRAPH2BIGWIG(ConvBase):
    """Converts :term:`BEDGRAPH` format to :term:`BIGWIG` format

    Conversion is based on bedGraph2BigWig tool. Note that an 
    argument --chrom-sizes is required.

    """
    _default_method = 'ucsc'

    def __init__(self, infile, outfile): #, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`BEDGRAPH` file.
        :param str outfile: (optional) output :term:`BIGWIG` file
        """
        super(BEDGRAPH2BIGWIG, self).__init__(infile, outfile)

    @requires("bedGraphToBigWig")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert bedgraph file in bigwig format using ucsc tool.
        https://genome.ucsc.edu/goldenpath/help/bigWig.html

        http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
        """
        chrom_sizes = kwargs.get("chrom_sizes", None)
        if chrom_sizes is None:
            raise ValueError("Must provide --chrom-sizes option")

        cmd = 'bedGraphToBigWig {infile}  {chrom_sizes} {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile, 
            chrom_sizes=chrom_sizes)
        self.execute(cmd)


    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--chrom-sizes",
            default=None,
            help="a two-column file/URL: <chromosome name> <size in bases>",
        )


