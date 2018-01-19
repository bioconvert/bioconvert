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
from bioconvert import ConvBase, extensions

import colorlog
logger = colorlog.getLogger(__name__)

__all__ = ["SAM2BAM"]


class SAM2BAM(ConvBase):
    """
    command used::
        samtools view -Sbh
    """

    def _method_samtools(self, *args, **kwargs):
        """
        Do the conversion  sorted :term`SAM` -> :term:'BAM`
        The result of the conversion is stored in the outputfile
        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        # -S means ignored (input format is auto-detected)
        # -b means output to BAM format
        # -h means include header in SAM output
        cmd = "samtools view -Sbh {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

