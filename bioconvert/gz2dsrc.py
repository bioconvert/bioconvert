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

__all__ = ["GZ2DSRC"]


class GZ2DSRC(ConvBase):
    """Convert :term:`GZ` file to :term:`BZ2` file

    Some description.

    """
    input_ext = [".gz"]
    output_ext = [".dsrc"]

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GZ file
        :param str outfile: output BZ2 filename

        """
        super(GZ2DSRC, self).__init__(infile, outfile, *args, **kargs)

    def _method_pigz_dscrc(self, *args, **kwargs):
        """some description"""
        # check integrity
        # cmd = "pigz -p{threads} --test {input}"
        # shell(cmd)

        # conversion
        cmd = "pigz -d -c -p {threads} {input} | dsrc c -s -t{threads} {output}"
        self.execute(cmd.format(
            threads=self.threads,
            input=self.infile,
            output=self.outfile))

        # integrity output
        #cmd = "pbzip2 {output} -p{threads} --test"
        #shell(cmd)

        # use self.infile, self.outfile

