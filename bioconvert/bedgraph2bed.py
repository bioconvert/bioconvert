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
"""Convert :term:`BEDGRAPH` file to :term:`BIGWIG` file."""
import os
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BEDGRAPH2BED"]


class BEDGRAPH2BED(ConvBase):
    """Converts a sequence alignment in :term:`BEDGRAPH` format to :term:`BED` format (3 columns)

    ::

        chr19   49302000    4930205    -1
        chr19   49302005    4930210    1
    
    becomes::

        chr19   4930201    -1
        chr19   4930202    -1
        chr19   4930203    -1
        chr19   4930204    -1
        chr19   4930205    -1
        chr19   4930206    1
        chr19   4930207    1
        chr19   4930208    1
        chr19   4930209    1
        chr19   4930210    1

    """
    _default_method = 'python'

    def __init__(self, infile, outfile): 
        """.. rubric:: constructor

        :param str infile: input :term:`BEDGRAPH` file.
        :param str outfile: output :term:`BED` file
        """
        super(BEDGRAPH2BED, self).__init__(infile, outfile)

    def _method_python(self, *args, **kwargs):
        """
        Convert bedgraph file in bed .
        """
        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
                for i, line in enumerate(fin.readlines()):
                    chrom, start, end, score = line.split()
                    assert start<end
                    for this in range(int(start), int(end)+1):
                        fout.write("{}\t{}\t{}\n".format(chrom, this, score))


