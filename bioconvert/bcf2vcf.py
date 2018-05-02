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
"""Convert :term:`BCF` file to :term:`VCF` file"""
from bioconvert import ConvBase

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)



class BCF2VCF(ConvBase):
    """Convert :term:`BCF` file to :term:`VCF` file

    """

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input BCF file
        :param str outfile: output VCF file

        command used::

            bcftools view 
        """
        super(BCF2VCF, self).__init__(infile, outfile, *args, **kargs)

    @requires("bcftools")
    def _method_bcftools(self, *args, **kwargs):

        # -O, --output-type b|u|z|v Output compressed BCF (b), uncompressed BCF
        # (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when
        # piping between bcftools subcommands to speed up performance
        cmd = "bcftools view {} -O v -o {}".format(self.infile, self.outfile)
        self.execute(cmd)



