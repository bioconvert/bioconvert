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
"""Convert :term:`VCF` file to :term:`BED` file"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


__all__ = ["VCF2BED"]


class VCF2BED(ConvBase):
    """
    Convert VCF file to BED file by extracting positions.

    Available methods:

    - awk::

        awk '! /\#/' INPUT | awk '{if(length($4) > length($5)) 
        print $1"\t"($2-1)"\t"($2+length($4)-1);
        else print $1"\t"($2-1)"\t"($2+length($5)-1)}' > OUTPUT

        This method report an interval of 1 for SNP, the length of the insertion or the length of 
        the deleted part in case of deletion.
    """
    _default_method = "awk"

    @requires("awk")
    def _method_awk(self, *args, **kwargs):
        """
        do the conversion :term`VCF` -> :term:'BED` using awk

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        awkcmd = """awk '{{if(length($4) > length($5)) print $1,($2-1),($2+length($4)-1); else print $1,($2-1),($2+length($5)-1)}}' OFS='\t'"""
        cmd = "awk '! /\#/' {} | {} > {}".format(self.infile, awkcmd, self.outfile)
        self.execute(cmd)
