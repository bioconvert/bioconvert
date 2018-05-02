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
"""Convert :term:`VCF` file to :term:`BCF` file"""
from bioconvert import ConvBase

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


__all__ = ["VCF2BCF"]


class VCF2BCF(ConvBase):
    """

    """

    @requires("bcftools")
    def _method_bcftools(self, *args, **kwargs):
        """

        command used::

            bcftools view -Sb

        :param args:
        :param kwargs:
        :return:
        """
        # -S means ignored (input format is VCF)
        # -b output BCF instead of VCF
        #cmd = "bcftools view -Sb {} >  {}".format(self.infile, self.outfile)

        # -O, --output-type b|u|z|v Output compressed BCF (b), uncompressed BCF
        # (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when
        # piping between bcftools subcommands to speed up performance
        cmd = "bcftools view {} -O b -o {}".format(self.infile, self.outfile)
        self.execute(cmd)



