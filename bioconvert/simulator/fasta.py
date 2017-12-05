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
"""Naive Fasta simulator for testing"""
import textwrap


class FastaSim():
    def __init__(self, outfile):
        self.outfile = outfile
        self.nreads = 1000000
        # do we want a wrap or unwrap version ?
        # For now, unwrap
        self.wrap = False
        self.wrap_length = 80
        self.read_length = 250

    def simulate(self):
        RL = self.read_length
        with open(self.outfile, "w") as fout:
            for i in range(self.nreads):
                fout.write(">identifier whatever it means but long enough\n")
                sequence = "ACGT" * (RL // 4) + "A" * (RL % 4) + "\n"
                if self.wrap is False:
                    fout.write(sequence)
                else:
                    wrapped = textwrap.wrap(sequence, self.wrap_length)
                    fout.write("\n".join(wrapped))
