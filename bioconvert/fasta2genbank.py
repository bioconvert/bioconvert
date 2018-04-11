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

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

__all__ = ["FASTA2GENBANK"]


class FASTA2GENBANK(ConvBase):
    """Convert :term:`FASTA` file to :term:`GENBANK` file"""

    # squizz works as well but keeps lower cases while 
    # biopython uses upper cases
    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input FASTA file
        :param str outfile: output GENBANK filename

        """
        super(FASTA2GENBANK, self).__init__(infile, outfile, *args, **kargs)

    @requires("squizz")
    def _method_squizz(self, *args, **kwargs):
        """Header is less informative than the one obtained with biopython"""
        cmd = "squizz {} -f fasta -c genbank > {} ".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        print("Using DNA alphabet for now")
        from Bio import SeqIO, Alphabet
        SeqIO.convert(self.infile, "fasta", self.outfile, "genbank",
            alphabet=Alphabet.generic_dna)
