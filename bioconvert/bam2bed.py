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
"""Convert :term:`BAM` format to :term:`BED` formats"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)

__all__ = ["BAM2BED"]


class BAM2BED(ConvBase):
    """Convert sorted :term:`BAM` file into :term:`BED` file 

    Available methods:

    - samtools::

        samtools depth -aa INPUT > OUTPUT

    - bedtools::

        bedtools genomecov -d -ibam INPUT > OUTPUT


    .. plot::

         from bioconvert.bam2bed import BAM2BED
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".bed") as fh:
             infile = bioconvert_data("test_measles.sorted.bam")
             convert = BAM2BED(infile, fh.name)
             convert.boxplot_benchmark()


    Note that this BED format is of the form::

        chr1    1   0
        chr1    2   0
        chr1    3   0
        chr1    4   0
        chr1    5   0

    that is contig name, position, coverage

    .. warning:: the BED file must be sorted. This can be achieved with
        bamtools.
    """
    _default_method = "samtools"

    def __init__(self, infile, outfile):
        """.. rubric:: Constructor

        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """
        do the conversion sorted :term:`BAM` -> :term:`BED` using samtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "samtools depth -aa {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires("bedtools")
    def _method_bedtools(self, *args, **kwargs):
        """
        do the conversion sorted :term`BAM` -> :term:'BED` using bedtools

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "bedtools genomecov -d -ibam {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)
