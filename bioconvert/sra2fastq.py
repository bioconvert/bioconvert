###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""Convert :term:`SRA` format to :term:`FASTQ` format"""

import os
import shutil
import subprocess
import tempfile

from bioconvert import ConvBase
from bioconvert.core.decorators import requires


class SRA2FASTQ(ConvBase):
    """Download FASTQ from SRA archive

    ::

        bioconvert sra2fastq ERR043367

    This may take some times since the files are downloaded from SRA website.

    """

    #: Default value
    _default_method = "fasterq_dump"

    # If test: will take only the first 10 reads from the sra file
    def __init__(self, infile, outfile, test=False):
        """.. rubric:: constructor


        https://edwards.flinders.edu.au/fastq-dump/

        library used: sra-toolkit
        """
        super().__init__(infile, outfile)
        self.test = test

    @requires("fastq-dump")
    def _method_fastq_dump(self, *args, **kwargs):
        """Uses Sratoolkit (fastq-dump) to convert a sra file to fastq

        `Fastq-dump documentation <https://edwards.flinders.edu.au/fastq-dump/>`_"""
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
            cmd = "fastq-dump {} {} --split-files -O {} {}".format(testcmd, compresscmd, tmpdir, infile)
            self.execute(cmd)
            cmd = "mv {}/{}_1.fastq{} {}_1.fastq{}".format(tmpdir, inname, gzext, outbasename, gzext)
            self.execute(cmd)
            cmd = "mv {}/{}_2.fastq{} {}_2.fastq{}".format(tmpdir, inname, gzext, outbasename, gzext)
            self.execute(cmd)
        else:
            cmd = "fastq-dump {} {} -O {} {}".format(testcmd, compresscmd, tmpdir, infile)
            self.execute(cmd)
            cmd = "mv {}/{}.fastq{} {}".format(tmpdir, inname, gzext, self.outfile)
            self.execute(cmd)
        shutil.rmtree(tmpdir)

    @requires("fasterq-dump")
    def _method_fasterq_dump(self, *args, **kwargs):
        """Uses fasterq-dump from SRA toolkit (modern replacement for fastq-dump).

        Handles both paired-end and single-end SRA files automatically by
        using ``--split-files`` and checking which output files were created.

        `fasterq-dump documentation <https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump>`_"""
        inname = os.path.split(os.path.splitext(self.infile)[0])[1]
        outbasename, ext = os.path.splitext(self.outfile)
        gzext = ""
        if ext == ".gz":
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

        # fasterq-dump with --split-files handles both paired and single-end;
        # for paired-end it creates <name>_1.fastq and <name>_2.fastq,
        # for single-end it creates <name>.fastq
        cmd = "fasterq-dump {} --split-files -O {} {}".format(testcmd, tmpdir, infile)
        self.execute(cmd)

        f1 = os.path.join(tmpdir, "{}_1.fastq".format(inname))
        f2 = os.path.join(tmpdir, "{}_2.fastq".format(inname))

        if os.path.exists(f1) and os.path.exists(f2):
            # Paired-end output
            if gzext:
                self.execute("gzip {}".format(f1))
                self.execute("gzip {}".format(f2))
            cmd = "mv {}_1.fastq{} {}_1.fastq{}".format(
                os.path.join(tmpdir, inname), gzext, outbasename, gzext
            )
            self.execute(cmd)
            cmd = "mv {}_2.fastq{} {}_2.fastq{}".format(
                os.path.join(tmpdir, inname), gzext, outbasename, gzext
            )
            self.execute(cmd)
        else:
            # Single-end output
            f = os.path.join(tmpdir, "{}.fastq".format(inname))
            if gzext:
                self.execute("gzip {}".format(f))
            cmd = "mv {}.fastq{} {}".format(
                os.path.join(tmpdir, inname), gzext, self.outfile
            )
            self.execute(cmd)
        shutil.rmtree(tmpdir)

    def isPairedSRA(self, filename):
        try:
            contents = subprocess.check_output(["fastq-dump", "-X", "1", "-Z", "--split-spot", filename])
        except subprocess.CalledProcessError:
            raise Exception("Error running fastq-dump on", filename)

        if contents.count(b"\n") == 4:
            return False
        elif contents.count(b"\n") == 8:
            return True
        else:
            raise Exception("Unexpected output from fast-dump on ", filename)
