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
"""Convert :term:`BAM` format to :term:`FASTQ` foarmat"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.utils import get_extension
import subprocess
import os


class BAM2FASTQ(ConvBase):
    """Convert sorted :term:`BAM` file into :term:`FASTQ` file

    Methods available are based on samtools [SAMTOOLS]_ or bedtools [BEDTOOLS]_.

    .. warning:: Using the bedtools method, the R1 and R2 reads must be next to 
        each other so that the reads are sorted similarly

    .. warning:: there is no guarantee that the R1/R2 output file are sorted
        similarly in paired-end case due to supp and second reads

    """
    _default_method = "samtools"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        """
        super(BAM2FASTQ, self).__init__(infile, outfile)

    """@requires("bamtools")
    def __method_bamtools(self, *args, **kwargs):

        # this method contains supplementary reads and we don't know 
        # what to do with them for now. So, this method is
        # commented. Indeed final R1 and R2 files will not be paired.

        cmd = "bamtools convert -format fastq -in {0} -out {1}".format(
            self.infile, self.outfile
        )
        self.execute(cmd)
    """

    @requires("bedtools")
    def _method_bedtools(self, *args, **kwargs):
        """Do the conversion :term:`BAM` -> :term:`Fastq` using bedtools

        """
        outbasename = os.path.splitext(self.outfile)[0]

        cmd = "bedtools bamtofastq -i {} -fq {}".format(self.infile, self.outfile)
        self.execute(cmd)

        output_ext = get_extension(self.outfile, remove_compression=True)

        # Due to the IO, paired reads are not always consecutive.
        # So, checking the first and second reads for paired data does not
        # work. Instead, we check the first 10 reads and check whether we have
        # at least one paired data. 

        data = []
        with open(self.outfile, "r") as fin:
            count = 0
            for i in range(40000):
                x = fin.readline()
                if len(x) == 0:
                    break
                data.append(x)

        from collections import Counter
        isPaired = 2 in Counter([x for i,x in enumerate(data) if i%4 == 0]).values()
        if isPaired:
            cmd = "bedtools bamtofastq -i {} -fq {}_1.{} -fq2 {}_2.{}".format(
                self.infile, outbasename, output_ext, outbasename, output_ext)
            self.execute(cmd)

            # Compress the output if required. We do not use compressor
            # since we may have two outputs.
            comp_ext = get_extension(self.outfile, remove_compression=False)
            if comp_ext in [".gz", ".dsrc", "bz2"]:
                from bioconvert.core.utils import compressor
                compressor("{}_1.{}".format(outbasename, output_ext), comp_ext)
                compressor("{}_2.{}".format(outbasename, output_ext), comp_ext)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """Do the conversion :term:`BAM` -> :term:`FASTQ` using samtools

        """
        cmd = "samtools fastq {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
        # Test if input bam file is paired
        p = subprocess.Popen("samtools view -c -f 1 {}".format(
            self.infile).split(),stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, universal_newlines=True)
        isPaired = p.communicate()[0].strip()

        # Collect the extension
        ext = os.path.splitext(self.outfile)[1]

        # FIXME: this compression code may be factorised ?

        output_ext = get_extension(self.outfile, remove_compression=True)

        # If the output file extension is compress extension
        if ext in [".gz",".bz2",".dsrc"]:
            outbasename = os.path.splitext(self.outfile)[0].split(".",1)[0]

            if ext == ".gz":
                compresscmd = "gzip -f"
            elif ext == ".bz2":
                compresscmd = "pbzip2 -f"
            else:
                compresscmd = "dsrc c"

            # When the input file is not paired and the output file needs to be compressed
            if isPaired == "0":
                cmd = "samtools fastq {} > {}.{}".format(self.infile,
                    outbasename, output_ext)
                self.execute(cmd)
                if ext == ".dsrc":
                    cmd = "{} {}.{} {}.{}.dsrc".format(compresscmd, outbasename,
                        output_ext, outbasename, output_ext)

                else:
                    cmd = "{} {}.{}".format(compresscmd, outbasename, output_ext)
                self.execute(cmd)
            # When the input file is paired and the output file needs to be compressed
            else:

                cmd = "samtools fastq -1 {}_1.{} -2 {}_2.{} -n {} ".format(outbasename, 
                    output_ext, outbasename, output_ext, self.infile)
                self.execute(cmd)
                if ext == ".dsrc":
                    cmd = "{} {}_1.{} {}_1.{}.dsrc".format(compresscmd,
                        outbasename, output_ext, outbasename, output_ext)
                    self.execute(cmd)
                    cmd = "{} {}_2.{} {}_2.{}.dsrc".format(compresscmd,
                        outbasename, output_ext, outbasename, output_ext)
                    self.execute(cmd)
                else:
                    cmd = "{} {}_1.{}".format(compresscmd, outbasename, output_ext)
                    self.execute(cmd)
                    cmd = "{} {}_2.{}".format(compresscmd, outbasename, output_ext)
                    self.execute(cmd)
        else:
            outbasename = os.path.splitext(self.outfile)[0]

            # When the input file is not paired
            if isPaired == "0":
                cmd = "samtools fastq {} > {}".format(self.infile, self.outfile)
                self.execute(cmd)
            # When the input file is paired
            else:
                #os.remove(self.outfile)
                cmd = "samtools fastq -1 {}_1.{} -2 {}_2.{} -n {} ".format(
                    outbasename, output_ext, outbasename, output_ext, self.infile)
                self.execute(cmd)
