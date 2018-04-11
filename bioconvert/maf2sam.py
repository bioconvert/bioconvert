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

import colorlog
logger = colorlog.getLogger(__name__)

__all__ = ["MAF2SAM"]


class MAF2SAM(ConvBase):
    """This is the Multiple alignment format or MIRA assembly format

    This is not Mutation Annotation Format (somatic)


    pbsim creates this kind of data

    Some references:

    - https://github.com/peterjc/maf2sam/
    - https://github.com/arq5x/nanopore-scripts/master/maf-convert.py
    - http://bioperl.org/formats/alignment_formats/MAF_multiple_alignment_format.html

    Those two codes were in Py2 at the time of this implementation. We re-used
    some of the information from maf-convert but the code in
    bioconvert.utils.maf can be considered original. 
    """

    def __init__(self, infile, outfile):
        super().__init__(infile, outfile)

    def _method_python(self, *args, **kwargs):
        from bioconvert.utils import maf
        conv = maf.MAF(self.infile, self.outfile)
        conv.to_sam()

