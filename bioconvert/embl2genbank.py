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
"""Convert :term:`EMBL` file to :term:`GENBANK` file"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

__all__ = ["EMBL2GENBANK"]


class EMBL2GENBANK(ConvBase):
    """Convert :term:`EMBL` file to :term:`GENBANK` file"""
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input EMBL file
        :param str outfile: output GENBANK filename

        """
        super(EMBL2GENBANK, self).__init__(infile, outfile, *args, **kargs)

    @requires(external_binary="squizz")
    def _method_squizz(self, *args, **kwargs):
        """Header is less informative than the one obtained with biopython"""
        cmd = "squizz {} -f embl -c genbank > {} ".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        from Bio import SeqIO
        SeqIO.convert(self.infile, "embl", self.outfile, "genbank")
