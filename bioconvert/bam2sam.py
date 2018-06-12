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
"""Convert :term:`SAM` file to :term:`BAM` file"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

import colorlog

logger = colorlog.getLogger(__name__)


class BAM2SAM(ConvBase):
    """

    .. plot::

        from bioconvert.bam2sam import BAM2SAM
        from bioconvert import bioconvert_data
        from easydev import TempFile

        with TempFile(suffix=".sam") as fh:
            infile = bioconvert_data("test_measles.sorted.bam")
            convert = BAM2SAM(infile, fh.name)
            convert.boxplot_benchmark()


    methods available:

    - samtools:
    - pysam
    - sambamba

    Could be implemented but not on bioconda:

    - sam-to-bam: Ogasawara T, Cheng Y, Tzeng T-HK (2016) Sam2bam:
        High-Performance Framework for NGS Data Preprocessing Tools. PLoS ONE
        11(11): e0167100. doi:10.1371/journal.pone.0167100

    """
    _default_method = "samtools"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        command used::

            samtools view -Sbh
        """
        super(BAM2SAM, self).__init__(infile, outfile, *args, **kargs)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        # -S means ignored (input format is auto-detected)
        # -h means include header in SAM output
        cmd = "samtools view -Sh {} -O SAM -o {}".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="pysam", external_binary="samtools")
    def _method_pysam(self, *args, **kwargs):
        import pysam
        pysam.sort("-o", self.outfile, self.infile)

    @requires("sambamba")
    def _method_sambamba(self, *args, **kwargs):
        cmd = "sambamba view {} -o {}".format(self.infile, self.outfile)
        self.execute(cmd)

