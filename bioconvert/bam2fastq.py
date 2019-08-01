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
"""Convert :term:`BAM` format to :term:`fastq` file"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires
import subprocess
import os


class BAM2FASTQ(ConvBase):
    """    Convert sorted :term:`BAM` file into :term:`FASTQ` file

    Methods available are based on samtools [SAMTOOLS]_ or bedtools [BEDTOOLS]_.

        .. warning:: Using the bedtools method, the R1 and R2 reads must be next each other
            so that the file to be properly separated

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

        # this method contains supplementary reads and we don't know what to do with them for now. So, this method is
        # commented. Indeed final R1 and R2 files will not be paired.

        cmd = "bamtools convert -format fastq -in {0} -out {1}".format(
            self.infile, self.outfile
        )
        self.execute(cmd)

    @requires("bedtools")
    def _method_bedtools(self, *args, **kwargs):
        """Do the conversion :term:`BAM` -> :term:`Fastq` using bedtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        outbasename = os.path.splitext(self.outfile)[0]

        cmd = "bedtools bamtofastq -i {} -fq {}".format(self.infile, self.outfile)
        self.execute(cmd)

        with open(self.outfile,'r') as paired_end:
            for i, line in enumerate(paired_end):
                if i == 0:
                    line_0 = line
                elif i == 4:
                    line_4 = line
                elif i > 4:
                    break
        if line_0 == line_4:
            cmd = "bedtools bamtofastq -i {} -fq {}_1.fastq -fq2 {}_2.fastq".format(self.infile, outbasename,outbasename)
            self.execute(cmd)
            os.remove(self.outfile)



    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """Do the conversion :term:`BAM` -> :term:`Fastq` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "samtools fastq {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
        # Test if input bam file is paired
        p = subprocess.Popen("samtools view -c -f 1 {}".format(
            self.infile).split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE,universal_newlines=True)
        isPaired =p.communicate()[0].strip()
        # Collect the extension
        ext = os.path.splitext(self.outfile)[1]

        # If the output file extension is compress extension
        if ext in [".gz",".bz2",".dsrc"]:
            outbasename = os.path.splitext(self.outfile)[0].split(".",1)[0]

            if ext == ".gz":
                compresscmd = "gzip"
            elif ext == ".bz2":
                compresscmd = "pbzip2 -f"
            else:
                compresscmd = "dsrc c"

            # When the input file is not paired and the output file needs to be compressed
            if isPaired == "0":
                cmd = "samtools fastq {} > {}.fastq".format(self.infile, outbasename)
                self.execute(cmd)
                if ext == ".dsrc":
                    cmd = "{} {}.fastq {}.fastq.dsrc".format(compresscmd, outbasename,outbasename)
                else:
                    cmd = "{} {}.fastq".format(compresscmd,outbasename)
                self.execute(cmd)
            # When the input file is paired and the output file needs to be compressed
            else:
                os.remove(self.outfile)
                cmd = "samtools fastq -1 {}_1.fastq -2 {}_2.fastq -n {} ".format(outbasename, outbasename, self.infile)
                self.execute(cmd)
                if ext == ".dsrc":
                    cmd = "{} {}_1.fastq {}_1.fastq.dsrc".format(compresscmd,outbasename, outbasename)
                    self.execute(cmd)
                    cmd = "{} {}_2.fastq {}_2.fastq.dsrc".format(compresscmd,outbasename, outbasename)
                    self.execute(cmd)
                else:
                    cmd = "{} {}_1.fastq".format(compresscmd,outbasename)
                    self.execute(cmd)
                    cmd = "{} {}_2.fastq".format(compresscmd,outbasename)
                    self.execute(cmd)


        else:
            outbasename = os.path.splitext(self.outfile)[0]

            # When the input file is not paired
            if isPaired == "0":
                cmd = "samtools fastq {} > {}".format(self.infile, self.outfile)
                self.execute(cmd)
            # When the input file is paired
            else:
                os.remove(self.outfile)
                cmd = "samtools fastq -1 {}_1.fastq -2 {}_2.fastq -n {} ".format(outbasename, outbasename, self.infile)
                self.execute(cmd)
