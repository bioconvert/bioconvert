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
"""ASQG validator"""
import colorlog
_log = colorlog.getLogger(__name__)


class ASQGLint():
    """ASQG format validator

    """
    def __init__(self, filename):
        self.filename = filename


    def validate(self):
        # read line by line. Checks 
        # - lines start with HT, VT or ED
        # - lines must be tab delimited
        # - lines VT must have 2 fields only
        with open(self.filename, "r") as fh:
            for i, line in enumerate(fh.readlines()):
                if line[0:2] == "HT":
                    pass
                elif line[0:2] == "VT":
                    if len(line.split()) < 3:
                        raise ValueError("Incorrect format on line %s. " % i +
                            "Line starting with VT must have at least 3 fields")
                elif line[0:2] == "ED":
                    if len(line.split()) != 11:
                        raise ValueError("Incorrect format on line %s. " % i +
                            "Line starting with ED must have 10 fields + ED prefix")
                elif len(line.strip()) == 0:
                    _log.warning("Found empty line on line %s" % line)
                else:
                    raise ValueError()

