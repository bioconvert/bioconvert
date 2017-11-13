# -*- coding: utf-8 -*-

##############################################################################
#  This file is part of Biokit software
#
#  Copyright (c) 2017 - Biokit Development Team
#
#  Distributed under the terms of the 3-clause BSD license.
#  The full license is in the LICENSE file, distributed with this software.
#
#  website: https://github.com/biokit/bioconvert
#  documentation: http://bioconvert.readthedocs.io
##############################################################################
""".. rubric:: Standalone application dedicated to conversion"""
import os
import sys
import colorlog
_log = colorlog.getLogger('bioconvert')

from bioconvert.core.registry import Registry
from bioconvert.core.utils import get_extension as getext


class Bioconvert(object):
    """Universal converter used by the standalone

    ::

        from bioconvert import Bioconvert
        c = Bioconvert("test.fastq", "test.fasta")


    """

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile: The path of the input file.
        :param str outfile: The path of The output file

        """
        if os.path.exists(infile) is False:
            msg = "Incorrect input file: %s" % infile
            _log.error(msg)
            raise ValueError(msg)

        self.inext = getext(infile)
        self.outext = getext(outfile)

        self.mapper = Registry()

        # From the input parameters 1 and 2, we get the module name
        try:
            _log.info("Input: {}".format(self.inext))
            _log.info("Output: {}".format(self.outext))
            class_converter = self.mapper[(self.inext, self.outext)]
            self.name = class_converter.__name__
        except KeyError:
            print(self.mapper)
            print(self.inext)
            print(self.outext)

            # Is the module name available in bioconvert ? If not, let us tell the user
            msg = "Requested input format ({}) to output format (({}) is not available in bioconvert"
            _log.critical(msg.format(self.inext, self.outext))
            _log.critical("Use --formats to know the available formats")
            sys.exit(1)

        self.converter = class_converter(infile, outfile)

    def __call__(self, *args, **kwargs):
        self.converter(*args, **kwargs)

    def boxplot_benchmark(self, *args, **kwargs):
        self.converter.boxplot_benchmark(*args,**kwargs)
