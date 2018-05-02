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
"""Convert :term:`BAM` format to :term:`fastq` file"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires


class BAM2Fastq(ConvBase):
    """Converts BAM 2 FastQ file

    .. warning:: the R1 and R2 reads are saved in the same file. Besides,
        there is no check that the read R1 and R2 alternates

    """
    _default_method = "bamtools"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor
        :param str infile:
        :param str outfile:
        library used::
            pysam (samtools)
        """
        super().__init__(infile, outfile)

    @requires("bamtools")
    def _method_bamtools(self, *args, **kwargs):
        # This fails with unknown error
        #pysam.bam2fq(self.infile, save_stdout=self.outfile)

        #cmd = "samtools fastq %s >%s" % (self.infile, self.outfile)
        #self.execute(cmd)

        # !!!!!!!!!!!!!!!!!! pysam.bam2fq, samtools fastq and bamtools convert
        # give differnt answers...

        cmd = "bamtools convert -format fastq -in {0} -out {1}".format(
            self.infile, self.outfile
        )
        self.execute(cmd)

    @requires("bedtools")
    def _method_bedtools(self, *args, **kwargs):
        """
        do the conversion :term`BAM` -> :term:'Fastq` using bedtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "bedtools bamtofastq -i {} -fq {}".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """
        do the conversion :term`BAM` -> :term:'Fastq` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "samtools fastq {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
