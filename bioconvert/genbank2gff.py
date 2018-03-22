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
""" Converts Genbank format to GFF3 using biocode
https://github.com/jorvis/biocode/

We may want to do it directly in python in the future,
without calling the external script
"""

from bioconvert import ConvBase

__all__ = ["GENBANK2GFF"]


class GENBANK2GFF(ConvBase):
    """Convert :term:`GENBANK` file to :term:`GFF` file

    Some description.

    """
    _default_method = "biocode"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GENBANK file
        :param str outfile: output GFF filename

        """
        super().__init__(infile, outfile)


    def _method_biocode(self, *args, **kwargs):
        """Uses scripts from biocode
        See: https://github.com/jorvis/biocode/
        https://github.com/jorvis/biocode/blob/master/gff/convert_genbank_to_gff3.py
        """
        cmd = "convert_genbank_to_gff3.py -i {} -o {} --no_fasta".format(
            self.infile, self.outfile)
        self.execute(cmd)
