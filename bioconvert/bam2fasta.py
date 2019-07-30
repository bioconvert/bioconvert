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
import subprocess
import os


class BAM2FASTA(ConvBase):
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
        do the conversion :term:`BAM` -> :term:`Fasta` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.

        .. note:: fasta are on one line
        """

        # Test if input bam file is paired
        p = subprocess.Popen("samtools view -c -f 1 {}".format(
            self.infile).split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
        isPaired =p.communicate()[0].strip()

        # Collect the extension
        ext = os.path.splitext(self.outfile)[1]

        # If the output file extension is compress extension
        if ext in [".gz",".bz2"]:
            outbasename = os.path.splitext(self.outfile)[0].split(".",1)[0]

            if ext == ".gz":
                compresscmd = "gzip"
            if ext == ".bz2":
                compresscmd = "pbzip2 -f"
            # When the input file is not paired and the output file needs to be compressed
            if isPaired == "0":
                cmd = "samtools fasta {} > {}.fasta".format(self.infile, outbasename)
                self.execute(cmd)
                cmd = "{} {}.fasta".format(compresscmd,outbasename)
                self.execute(cmd)
            # When the input file is paired and the output file needs to be compressed
            else:
                cmd = "samtools fasta -1 {}_1.fasta -2 {}_2.fasta -n {} ".format(outbasename, outbasename, self.infile)
                self.execute(cmd)
                cmd = "{} {}_1.fasta".format(compresscmd,outbasename)
                self.execute(cmd)
                cmd = "{} {}_2.fasta".format(compresscmd,outbasename)
                self.execute(cmd)

        else:
            outbasename = os.path.splitext(self.outfile)[0]

            # When the input file is not paired
            if isPaired == "0":
                cmd = "samtools fasta {} > {}".format(self.infile, self.outfile)
                self.execute(cmd)
            # When the input file is paired
            else:
                cmd = "samtools fasta -1 {}_1.fasta -2 {}_2.fasta -n {} ".format(outbasename, outbasename, self.infile)
                self.execute(cmd)

