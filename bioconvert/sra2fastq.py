# -*- coding: utf-8 -*-
#
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
#
##############################################################################
"""Convert :term:`SRA` format to :term:`Fastq` file"""
from bioconvert import ConvBase, extensions
import subprocess
import os

class Sra2Fastq(ConvBase):
    """Converts Sra 2 Fastq file

    """

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        library used::
            sra-toolkit
        """
        super().__init__(infile, outfile)
        self._default_method = "sratoolkit"

    def _method_sratoolkit(self, *args, **kwargs):
        """
        Uses Sratoolkit (fastq-dump) to convert a sra file to fastq
        """
        outname='{}.{}'.format(os.path.splitext(self.infile)[0],"fastq")
        inbasename=os.path.splitext(self.infile)[0]
        outbasename=os.path.splitext(self.outfile)[0]
        infile = self.infile

        # If the file does not exist locally, we take the basename
        # it should correspond to a SRA ID
        if os.path.isfile(infile) is False:
            infile = os.path.splitext(self.infile)[0]
            infile = os.path.split(infile)[1]

        if self.isPairedSRA(infile):
           cmd = "fastq-dump --split-files "+infile
           self.execute(cmd)
           if self.outfile!=outname:
                   cmd = "mv {0} {1}".format(inbasename+"_1.fastq", outbasename+"_1.fastq")
                   self.execute(cmd)
                   cmd = "mv {0} {1}".format(inbasename+"_2.fastq", outbasename+"_2.fastq")
                   self.execute(cmd)
        else:
           cmd = "fastq-dump {}".format(infile)
           self.execute(cmd)
           if self.outfile!=outname :
                cmd = "mv {0} {1}".format(inbasename+".fastq", self.outfile)
                self.execute(cmd)


    def isPairedSRA(self,filename):
        try:
            contents = subprocess.check_output(["fastq-dump","-X","1","-Z","--split-spot", filename]);
        except subprocess.CalledProcessError:
            raise Exception("Error running fastq-dump on",filename);

        if(contents.count(b'\n') == 4):
            return False;
        elif(contents.count(b'\n') == 8):
            return True;
        else:
            raise Exception("Unexpected output from fast-dump on ", filename);
        
