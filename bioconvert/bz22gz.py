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
from bioconvert import ConvBase, extensions

import colorlog
logger = colorlog.getLogger(__name__)


__all__ = ["BZ22GZ"]


class BZ22GZ(ConvBase):
    """Convert :term:`BZ2` file to :term:`GZ` file

    Some description.

    """

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input BZ2 file
        :param str outfile: output GZ filename

        """
        super(BZ22GZ, self).__init__(infile, outfile, *args, **kargs)

        self._default_method = "bz2_gz"

    def _method_bz2_gz(self, *args, **kwargs):
        # conversion
        cmd = "bunzip2 -c {input} | gzip > {output}".format(
            input=self.infile,
            output=self.outfile)
        self.execute(cmd)
