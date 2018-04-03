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
"""BEDGRAPH2BIGWIG conversion """
import os
import colorlog

from bioconvert.core.base import ConvArg, ConvBase

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BEDGRAPH2BIGWIG"]


class BEDGRAPH2BIGWIG(ConvBase):
    """Converts a sequence alignment in :term:`BEDGRAPH` format to :term:`BIGWIG` format

    Conversion is based on ucsc bedGraph2BigWig tool

    """
    _default_method = 'ucsc'

    def __init__(self, infile, outfile): #, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`BEDGRAPH` file.
        :param str outfile: (optional) output :term:`BIGWIG` file
        """
        super(BEDGRAPH2BIGWIG, self).__init__(infile, outfile)

    @requires("bedGraphToBigWig")
    def _method_ucsc(self, *args, **kwargs):
        """
        Convert bedgraph file in bigwig format using ucsc tool.
        https://genome.ucsc.edu/goldenpath/help/bigWig.html
        
        http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
        """
        chrom_sizes = kwargs.get("chrom_sizes", None)
        if chrom_sizes is None:
            raise ValueError("Must provide --chrom-sizes option")

        cmd = 'bedGraphToBigWig {infile}  {chrom_sizes} {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile, 
            chrom_sizes=chrom_sizes)
        self.execute(cmd)


    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--chrom-sizes",
            default=None,
            help="a two-column file/URL: <chromosome name> <size in bases>",
        )


