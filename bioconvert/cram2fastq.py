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
"""Convert :term:`CRAM` file to :term:`FASTQ` format"""
from bioconvert import ConvBase
import os
from bioconvert.core.utils import get_extension
from easydev.multicore import cpu_count
import subprocess

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


class CRAM2FASTQ(ConvBase):
    """Convert :term:`CRAM` file to :term:`FASTQ` file

    Methods available are based on samtools [SAMTOOLS]_.

    """
    _default_method = "samtools"
    _threading = True

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input CRAM file
        :param str outfile: output FASTQ filename

        """
        super(CRAM2FASTQ, self).__init__(infile, outfile, *args, **kargs)

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
                cmd = "samtools fastq -@ {} {} > {}.{}".format(self.threads, 
                    self.infile, outbasename, output_ext)
                self.execute(cmd)
                if ext == ".dsrc":
                    cmd = "{} {}.{} {}.{}.dsrc".format(compresscmd, 
                        outbasename, output_ext, outbasename, output_ext)
                else:
                    cmd = "{} {}.{}".format(compresscmd, outbasename, output_ext)
                self.execute(cmd)
            # When the input file is paired and the output file needs to be compressed
            else:

                cmd = "samtools fastq -@ {} -1 {}_1.{} -2 {}_2.{} -n {} ".format(
                    self.threads, outbasename, output_ext, outbasename, output_ext, 
                    self.infile)
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
                cmd = "samtools fastq -@ {} {} > {}".format(self.threads, 
                    self.infile, self.outfile)
                self.execute(cmd)
            # When the input file is paired
            else:
                cmd = "samtools fastq -@ {} -1 {}_1.{} -2 {}_2.{} -n {} ".format(self.threads, 
                    outbasename, output_ext, outbasename, output_ext, self.infile)
                self.execute(cmd)
