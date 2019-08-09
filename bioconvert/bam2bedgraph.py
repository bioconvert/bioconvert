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
"""Convert :term:`BAM` format to :term:`BEDGRAPH` format"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BAM2BEDGRAPH"]


class BAM2BEDGRAPH(ConvBase):
    """Convert sorted :term:`BAM` file into :term:`BEDGRAPH` file

    Compute the coverage (depth) in BEDGRAPH.
    Regions with zero coverage are also reported.


    Note that this BEDGRAPH format is of the form::

        chrom chromStart chromEnd dataValue

    Note that consecutive positions with same values are compressed.

    ::

        chr1    0   75  0
        chr1    75  176 1
        chr1    176  177 2


    .. warning:: the BAM file must be sorted. This can be achieved with
        bamtools.


    Methods available are based on bedtools [BEDTOOLS]_ and mosdepth
    [MOSDEPTH]_.
    """
    # 4 minutes with bedtools and 20s with mosdepth
    _default_method = "bedtools"
    _threading = True


    def __init__(self, infile, outfile):
        """.. rubric:: Constructor

        :param str infile: The path to the input BAM file.
            **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)

    @requires("bedtools")
    def _method_bedtools(self, *args, **kwargs):
        """Do the conversion using bedtools"""
        cmd = "bedtools genomecov -bga -ibam {} > {}".format(self.infile,
                                                             self.outfile)
        self.execute(cmd)

    @requires("mosdepth")
    def _method_mosdepth(self, *args, **kwargs):
        """Do the conversion using mosdepth"""
        # For testing, we need to save into a specific temporary directory
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                cmd = "mosdepth {}/.bioconvert -t {}  {}".format(tmpdir, self.threads, self.infile)
                self.execute(cmd)

                if self.outfile.endswith(".gz"):
                    pass
                else:
                    cmd = "gunzip -c {}/.bioconvert.per-base.bed.gz > {}".format(tmpdir, self.outfile)
                    self.execute(cmd)
            except Exception as err:
                raise(err)
            finally:
                cmd = "rm -f {name}/.bioconvert.per-base.bed.gz {name}/.bioconvert.per-base.bed.gz.csi"
                cmd += " {name}/.bioconvert.mosdepth.global.dist.txt"
                self.execute(cmd.format(name=tmpdir))
