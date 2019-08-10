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
"""Convert :term:`ABI` format to :term:`QUAL` format"""
from bioconvert import ConvBase
from bioconvert import requires


__all__ = ["ABI2QUAL"]


class ABI2QUAL(ConvBase):
    """Convert :term:`ABI` file to :term:`QUAL` file

    :term:`ABI` files are created by ABI sequencing machine and
    includes PHRED quality scores for base calls. This allows
    the creation of :term:`QUAL` files.

    Method implemented is based on BioPython [BIOPYTHON]_.

    """

    _default_method = "biopython"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input ABI file
        :param str outfile: output QUAL filename

        """
        super(ABI2QUAL, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="biopython")
    def _method_biopython(self, *args, **kwargs):
        from Bio import SeqIO
        records = SeqIO.parse(self.infile, "abi")
        # output using SeqIO.write(records, self.outfile, "qual") is not
        # standard so we write our own conversion here below
        with open(self.outfile, "w") as fout:
            for rec in records:
                header = rec.name
                qual = rec.letter_annotations['phred_quality']
                qual = "".join([str(x) for x in qual])
                fout.write(">{}\n".format(header))
                fout.write("{}\n".format(qual))
