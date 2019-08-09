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
"""Convert :term:`SRA` format to :term:`FASTA` format"""

from bioconvert import ConvBase
import subprocess
import os
import tempfile
import shutil

from bioconvert.core.decorators import requires


class SRA2FASTQ(ConvBase):
    """Download FASTQ from SRA archive

    ::

        bioconvert sra2fastq ERR043367

    This may take some times since the files are downloaded from SRA website.

    """
    _default_method = "sratoolkit"

    # If test: will take only the first 10 reads from the sra file
    def __init__(self, infile, outfile, test=False):
        """.. rubric:: constructor
    
        :param str infile:
        :param str outfile:

        library used: sra-toolkit
        """
        super().__init__(infile, outfile)
        self.test = test

    @requires("fastq-dump")
    def _method_sratoolkit(self, *args, **kwargs):
        """
        Uses Sratoolkit (fastq-dump) to convert a sra file to fastq
        """
        inname = os.path.split(os.path.splitext(self.infile)[0])[1]
        outbasename, ext = os.path.splitext(self.outfile)
        compresscmd = ""
        gzext = ""
        if ext == ".gz":
            compresscmd = "--gzip"
            gzext = ".gz"
            outbasename = os.path.splitext(outbasename)[0]

        infile = self.infile
        # If the file does not exist locally, we take the basename
        # it should correspond to a SRA ID
        if os.path.isfile(infile) is False:
            infile = inname

        tmpdir = tempfile.mkdtemp()
        testcmd = ""
        # If in test mode, we retrieve only 10 reads from sra
        if self.test:
            testcmd = "-X 10"
        if self.isPairedSRA(infile):
            cmd = "fastq-dump {} {} --split-files -O {} {}".format(
                testcmd, compresscmd, tmpdir, infile)
            self.execute(cmd)
            cmd = "mv {}/{}_1.fastq{} {}_1.fastq{}".format(
                tmpdir, inname, gzext, outbasename, gzext)
            self.execute(cmd)
            cmd = "mv {}/{}_2.fastq{} {}_2.fastq{}".format(
                tmpdir, inname, gzext, outbasename, gzext)
            self.execute(cmd)
        else:
            cmd = "fastq-dump {} {} -O {} {}".format(
                testcmd, compresscmd, tmpdir, infile)
            self.execute(cmd)
            cmd = "mv {}/{}.fastq{} {}".format(tmpdir,
                                               inname, gzext, self.outfile)
            self.execute(cmd)
        shutil.rmtree(tmpdir)

    def isPairedSRA(self, filename):
        try:
            contents = subprocess.check_output(
                ["fastq-dump", "-X", "1", "-Z", "--split-spot", filename])
        except subprocess.CalledProcessError:
            raise Exception("Error running fastq-dump on", filename)

        if(contents.count(b'\n') == 4):
            return False
        elif(contents.count(b'\n') == 8):
            return True
        else:
            raise Exception("Unexpected output from fast-dump on ", filename)
