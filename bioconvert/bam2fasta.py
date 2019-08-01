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
import itertools


class BAM2FASTA(ConvBase):
    """Bam2Fasta converter

    Convert sorted :term:`BAM` file into :term:`FASTA` file

    Methods available are based on samtools [SAMTOOLS]_.
    """
    _default_method = "samtools"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        library used: pysam (samtools)

        """
        super().__init__(infile, outfile)

    @requires("bamtools")
    def __method_bamtools(self, *args, **kwargs):
        """

        .. note:: fastq are split on several lines (80 characters)

        """
        # this method contains supplementary reads and we don't know what to do with them for now. So, this method is
        # commented. Indeed final R1 and R2 files will not be paired.

        cmd = "bamtools stats -in '%s' | sed '12!d' | awk '{print $3}' " % (self.infile)
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,universal_newlines=True)
        isPaired = ps.communicate()[0].strip()

        # Collect the extension
        ext = os.path.splitext(self.outfile)[1]

        cmd = "bamtools convert -format fasta -in {} -out {}".format(
            self.infile, self.outfile)
        self.execute(cmd)

        if isPaired != "0":
            with open(self.outfile, "r") as paired_end, open("Reads_1.fasta","w") as out1, open("Reads_2.fasta","w") as out2:
                outfile_cycler = itertools.cycle((out1, out2))
                for line in paired_end:
                    outfile = next(outfile_cycler)
                    outfile.write(line)
                    outfile.write(next(paired_end))


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

