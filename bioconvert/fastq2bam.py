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

"""Convert :term:`FASTQ` to :term:`BAM`"""
from bioconvert import ConvBase

from bioconvert.core.decorators import in_gz, requires


class FASTQ2BAM(ConvBase):
    """Convert :term:`FASTQ` to :term:`BAM`"""
    _default_method = "fastqutils"

    # infile: read 1, infile2: read2 if paired-end
    def __init__(self, infile, outfile, infile2=None, *args, **kwargs):
        """
        :param str infile: The path to the input FASTA file.
        :param str outfile: The path to the output file.
        """
        super().__init__(infile, outfile)
        # use readfq for now because pure python are fast enough
        # for production, could use seqtk which seems the fastest method though
        # Make sure that the default handles also the compresssion
        self.infile2 = infile2

    @requires("conda")
    @in_gz
    def _method_fastqutils(self, *args, **kwargs):
        """
        Converts a fastq file to an unaligned bam file
        """
        self.install_tool('fastqutils')
        if self.infile2 is not None:
            cmd = "fastqutils tobam -1 {} -2 {} -o {}".format(
                self.infile, self.infile2, self.outfile)
        else:
            cmd = "fastqutils tobam -1 {} -o {}".format(
                self.infile, self.outfile)
        self.execute(cmd)


#TODO: could use picard as follows:
# picard FastqToSam FASTQ=sd_0001.fastq OUTPUT=test2.bam  READ_GROUP_NAME=test
#     SAMPLE_NAME=test LIBRARY_NAME=sim PLATFORM=pacbio

