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
"""Convert :term:`CRAM` file to :term:`FASTQ` file"""
import os
from bioconvert import ConvBase
from easydev.multicore import cpu_count

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


class CRAM2FASTQ(ConvBase):
    """Convert :term:`CRAM` file to :term:`FASTQ` file

    The conversion requires the reference corresponding to the input file
    It can be provided as an argument in the constructor. Otherwise,
    a local file with same name as the input file but an .fa extension is
    looked for. Otherwise, we ask for the user to provide the input file.
    This is useful for the standalone application.

    """
    _default_method = "samtools"

    def __init__(self, infile, outfile, reference=None, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input FASTQ file
        :param str outfile: output filename
        :param str reference: reference file in :term:`FASTA` format

        command used::

            samtools view -@ <thread> -Sh -T <reference> in.cram > out.sam

        .. note:: the API related to the third argument may change in the future.
        """
        super(CRAM2FASTQ, self).__init__(infile, outfile, *args, **kargs)
        self.threads = cpu_count()

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        # -h means include header in FASTQ output
        # TODO: ideally we need to first use
        # samtools collate in.cram out.bam and then
        # samtools fastq -1 1.fastq -2 2.fastq out.bam
        cmd = "samtools fastq -@ {} -1 {} -2 {} {}".format(self.threads,
                self.outfile+".1", self.outfile+".2", self.infile)
        self.execute(cmd)
