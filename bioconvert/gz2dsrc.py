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
#
""" Convert a compressed fastq.gz file to :term:`DSRC` compression format """

from bioconvert import ConvBase

__all__ = ["GZ2DSRC"]


class GZ2DSRC(ConvBase):
    """Convert compressed fastq.gz file into `DSRC` compressed file

    .. plot::

         from bioconvert.gz2dsrc import GZ2DSRC
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".dsrc") as fh:
             infile = bioconvert_data("test_SP1.fq.gz")
             convert = GZ2DSRC(infile, fh.name)
             convert.boxplot_benchmark()

    """
    input_ext = [".gz"]
    output_ext = [".dsrc"]

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GZ filename
        :param str outfile: output DSRC filename

        """
        super(GZ2DSRC, self).__init__(infile, outfile, *args, **kargs)
        self._default_method = "pigzdsrc"

    def _method_pigzdsrc(self, *args, **kwargs):
        """
        do the conversion gz -> :term:'DSRC`

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "pigz -d -c -p {threads} {input} | dsrc c -s -t{threads} {output}"
        self.execute(cmd.format(
            threads=self.threads,
            input=self.infile,
            output=self.outfile))


