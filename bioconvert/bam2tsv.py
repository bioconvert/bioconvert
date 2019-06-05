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
"""Convert :term:`bam` file to :term:`tsv` file"""
import os
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

import colorlog

logger = colorlog.getLogger(__name__)

class BAM2TSV(ConvBase):
    """

    .. plot::

        from bioconvert.bam2tsv import BAM2TSV
        from bioconvert import bioconvert_data
        from easydev import TempFile

        with TempFile(suffix=".tsv") as fh:
            infile = bioconvert_data("test_measles.sorted.bam")
            convert = BAM2TSV(infile, fh.name)
            convert.boxplot_benchmark()


    methods available:

    - samtools
    - pysam
    """
    _default_method = "samtools"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        command used::

            samtools index && samtools idxstats
        """
        super(BAM2TSV, self).__init__(infile, outfile, *args, **kargs)


    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        with open(self.outfile, 'wt') as out:
            out.write("Reference sequence name\tSequence length\t"
                      "Mapped reads\tUnmapped reads{}".format(os.linesep))
        cmd = ("samtools index {0} && samtools idxstats {0} >> {1}"
                .format(self.infile, self.outfile))
        self.execute(cmd)

    @requires(python_library="pysam", external_binary="samtools")
    def _method_pysam(self, *args, **kwargs):
        import pysam
        # index the bam file
        pysam.index(self.infile)
        # create count table
        with open(self.outfile, 'wt') as out:
            out.write("Reference sequence name\tSequence length\t"
                      "Mapped reads\tUnmapped reads{}".format(os.linesep))
            for line in pysam.idxstats(self.infile):
                out.write(line) 
