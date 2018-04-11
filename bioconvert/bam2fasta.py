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
"""Convert :term:`BAM` format to :term:`FASTA` file"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires


class BAM2Fasta(ConvBase):
    """Bam2Fasta converter

    Wrapper of bamtools to convert bam file to fasta file.

    """
    _default_method = "bamtools"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        library used: pysam (samtools)
        """
        super().__init__(infile, outfile)

    @requires("bamtools")
    def _method_bamtools(self, *args, **kwargs):
        """

        .. note:: fastq are split on several lines (80 characters)

        """
        # Another idea is to use pysam.bam2fq but it fails with unknown error
        #pysam.bam2fq(self.infile, save_stdout=self.outfile)
        #cmd = "samtools fastq %s >%s" % (self.infile, self.outfile)
        #self.execute(cmd)
        # !!!!!!!!!!!!!!!!!! pysam.bam2fq, samtools fastq and bamtools convert
        # give differnt answers...

        cmd = "bamtools convert -format fasta -in {0} -out {1}".format(
            self.infile, self.outfile
        )
        self.execute(cmd)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """
        do the conversion :term`BAM` -> :term:'Fasta` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.

        .. note:: fasta are on one line
        """
        cmd = "samtools fasta {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

