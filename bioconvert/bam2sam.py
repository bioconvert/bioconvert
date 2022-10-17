###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""Convert :term:`SAM` file to :term:`BAM` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

import colorlog

logger = colorlog.getLogger(__name__)

__all__ = ["BAM2SAM"]


class BAM2SAM(ConvBase):
    """Convert :term:`BAM` file to :term:`SAM` file


    Methods available are based on samtools [SAMTOOLS]_ , sam-to-bam [SAMTOBAM]_ ,
    sambamba [SAMBAMBA]_ and pysam [PYSAM]_.

    """

    #: default value
    _default_method = "sambamba"
    _threading = True

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        """
        super(BAM2SAM, self).__init__(infile, outfile, *args, **kargs)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        # -S means ignored (input format is auto-detected)
        # -h means include header in SAM output
        """Here we use the SAMtools tool.

        `SAMtools documentation <http://www.htslib.org/doc/samtools.html>`_"""
        cmd = "samtools view -Sh {} --threads {} -O SAM -o {}"
        cmd = cmd.format(self.infile, self.threads, self.outfile)
        self.execute(cmd)

    @requires(python_library="pysam", external_binary="samtools")
    def _method_pysam(self, *args, **kwargs):
        """We use here the python module Pysam.

        `Pysam documentation <https://pysam.readthedocs.io/en/latest/api.html>`_"""
        import pysam

        pysam.sort("-o", self.outfile, self.infile)

    @requires("sambamba")
    def _method_sambamba(self, *args, **kwargs):
        """Here we use the Sambamba tool.
        This is the default method because it is the fastest.

        `Sambamba documentation <https://lomereiter.github.io/sambamba/docs/sambamba-view.html>`_"""
        # do not use --header that returns only the header...
        cmd = "sambamba view  {} -o {} -t {}"
        cmd = cmd.format(self.infile, self.outfile, self.threads)
        self.execute(cmd)
