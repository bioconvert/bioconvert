# -*- coding: utf-8 -*-
##############################################################################
#  This file is part of Bioconvert software
#
#  Copyright (c) 2017 - Bioconvert Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
""" description """
import colorlog
from bioconvert import ConvBase, extensions
# This has side effects in the registry !! 
#from bioconvert.bam2bedgraph import BAM2BEDGRAPH
#from bioconvert.bedgraph2bigwig import BEDGRAPH2BIGWIG
from easydev import TempFile

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)

__all__ = ["BAM2BIGWIG"]


class BAM2BIGWIG(ConvBase):
    """Convert :term:`BAM` file to :term:`BIGWIG` file

    Some description.

    """
    _default_method = "bamCoverage"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input BAM file
        :param str outfile: output BIGWIG filename

        command used::
            bamCoverage -bam {} –-outFileFormat bigwig --outFileName {}


        """
        super(BAM2BIGWIG, self).__init__(infile, outfile, *args, **kargs)



    @requires("bamCoverage")
    def _method_bamCoverage(self, *args, **kwargs):
        """run bam2bigwig from deeptools package"""
        cmd = "bamCoverage --bam {} --outFileFormat bigwig --outFileName {}".format(
                self.infile, self.outfile)
        self.execute(cmd)

    @requires("bedtools")
    def _method_ucsc(self, *args, **kwargs):
        """run bam2bigwig using bioconvert/bedtools bam2bedgraph and ucsc tool bedGraphToBigWig"""
        from bioconvert.bam2bedgraph import BAM2BEDGRAPH
        from bioconvert.bedgraph2bigwig import BEDGRAPH2BIGWIG
        with TempFile(suffix='.bedgraph') as fh:
            convertbam2bed = BAM2BEDGRAPH(self.infile, fh.name)
            convertbam2bed()
            convertbed2bw = BEDGRAPH2BIGWIG(fh.name, self.outfile)
            convertbed2bw()


