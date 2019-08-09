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

"""Convert :term:`BAM` file to :term:`BIGWIG` format"""

import colorlog
from bioconvert import ConvBase 
from bioconvert.core.base import ConvArg


from easydev import TempFile

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)

__all__ = ["BAM2BIGWIG"]


class BAM2BIGWIG(ConvBase):
    """Convert :term:`BAM` file to :term:`BIGWIG` file

    Convert BAM into a binary version of :term:`WIG` format.

    Methods are base on bamCoverage [DEEPTOOLS]_ and bedGraphToBigWig from
    wiggletools [WIGGLETOOLS]_. Wiggletools method requires an extra argument
    (--chrom-sizes) therefore default one is bamCoverage for now.

    Moreover, the two methods do not return exactly the same info!

    You can check this by using bioconvert to convert into a human readable file
    such as wiggle. We will use the bamCoverage as our default conversion.

    """
    _default_method = "bamCoverage"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input BAM file
        :param str outfile: output BIGWIG filename

        """
        super(BAM2BIGWIG, self).__init__(infile, outfile, *args, **kargs)

    @requires("bamCoverage")
    def _method_bamCoverage(self, *args, **kwargs):
        """run bamCoverage package"""
        cmd = "bamCoverage --bam {} --outFileFormat bigwig --outFileName {}".format(
                self.infile, self.outfile)
        self.execute(cmd)

    @requires(external_binaries=["bedGraphToBigWig", "bedtools"])
    def _method_ucsc(self, *args, **kwargs):
        """Run ucsc tool bedGraphToBigWig. 

        Requires extra argument (chrom_sizes) required by the bioconvert
        stanalone. """
        from bioconvert.bam2bedgraph import BAM2BEDGRAPH
        from bioconvert.bedgraph2bigwig import BEDGRAPH2BIGWIG

        chrom_sizes = kwargs.get("chrom_sizes", None)

        with TempFile(suffix='.bedgraph') as fh:
            convertbam2bed = BAM2BEDGRAPH(self.infile, fh.name)
            convertbam2bed()
            convertbed2bw = BEDGRAPH2BIGWIG(fh.name, self.outfile)
            convertbed2bw(chrom_sizes=chrom_sizes)

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--chrom-sizes",
            default=None,
            help="a two-column file/URL: <chromosome name> <size in bases>. "
                 "Used by the ucsc method only",
        )
