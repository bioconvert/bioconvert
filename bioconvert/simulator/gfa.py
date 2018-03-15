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
"""Naive GFA simulator for testing"""


class GFASim():
    """Simple GFA simulator




    """
    def __init__(self, outfile):
        self.outfile = outfile
        self.nreads = 1000000
        self.read_length = 250

    def simulate(self):
        RL = self.read_length
        with open(self.outfile, "w") as fout:
            for i in range(self.nreads):
                sequence = "ACGT" * (RL // 4) + "A" * (RL % 4) 
                fout.write("S contig{} {}\n".format(i, sequence))



