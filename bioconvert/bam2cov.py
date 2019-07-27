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
"""Convert :term:`BAM` format to :term:`COV` format"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)

__all__ = ["BAM2COV"]


class BAM2COV(ConvBase):
    """Convert sorted :term:`BAM` file into :term:`COV` file 


    Note that the COV format is of the form::

        chr1    1   0
        chr1    2   0
        chr1    3   0
        chr1    4   0
        chr1    5   0

    that is contig name, position, coverage.

    .. warning:: the BAM file must be sorted. This can be achieved with
        bamtools using *bamtools sort -in INPUT.bam*

    Methods available are based on samtools [SAMTOOLS]_ or bedtools [BEDTOOLS]_.
    """
    _default_method = "samtools"

    def __init__(self, infile, outfile):
        """.. rubric:: Constructor

        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """Do the conversion sorted :term:`BAM` -> :term:`BED` using samtools"""
        cmd = "samtools depth -aa {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires("bedtools")
    def _method_bedtools(self, *args, **kwargs):
        """Do the conversion sorted :term:`BAM` -> :term:`BED` using bedtools"""
        cmd = "bedtools genomecov -d -ibam {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
