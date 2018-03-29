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
"""Convert :term:`SAM` file to :term:`BAM` file"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

import colorlog

logger = colorlog.getLogger(__name__)


class BAM2SAM(ConvBase):
    """

    .. plot::

        from bioconvert.bam2sam import BAM2SAM
        from bioconvert import bioconvert_data
        from easydev import TempFile

        with TempFile(suffix=".sam") as fh:
            infile = bioconvert_data("test_measles.sorted.bam")
            convert = BAM2SAM(infile, fh.name)
            convert.boxplot_benchmark()


    methods available:

    - samtools:
    - pysam
    - sambamba

    Could be implemented but not on bioconda:

    - sam-to-bam: Ogasawara T, Cheng Y, Tzeng T-HK (2016) Sam2bam:
        High-Performance Framework for NGS Data Preprocessing Tools. PLoS ONE
        11(11): e0167100. doi:10.1371/journal.pone.0167100

    """
    _default_method = "samtools"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        command used::

            samtools view -Sbh
        """
        super(BAM2SAM, self).__init__(infile, outfile, *args, **kargs)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        # -S means ignored (input format is auto-detected)
        # -h means include header in SAM output
        cmd = "samtools view -Sh {} -O SAM -o {}".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="pysam", external_binary="samtools")
    def _method_pysam(self, *args, **kwargs):
        import pysam
        pysam.sort("-o", self.outfile, self.infile)

    @requires("sambamba")
    def _method_sambamba(self, *args, **kwargs):
        cmd = "sambamba view {} -o {}".format(self.infile, self.outfile)
        self.execute(cmd)

