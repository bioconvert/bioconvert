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
"""Convert :term:`BIGBED` format to :term:`WIGGLE` format"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BIGBED2WIGGLE"]


class BIGBED2WIGGLE(ConvBase):
    """Convert sorted :term:`BIGBED` file into :term:`WIGGLE` file

    Methods available are based on wiggletools [WIGGLETOOLS]_.
    """
    _default_method = "wiggletools"

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BIGBED file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super(BIGBED2WIGGLE, self).__init__(infile, outfile)

    @requires("wiggletools")
    def _method_wiggletools(self, *args, **kwargs):
        """

        """
        import os
        from easydev import TempFile

        # with need a unique name, that does not exists for the symlink 
        # Fixes #233
        fname = None
        with TempFile(suffix=".bb") as ftemp:
            fname = ftemp.name

        os.symlink(os.path.abspath(self.infile), ftemp.name)

        try:
            cmd = "wiggletools {} > {}".format(ftemp.name, self.outfile)
            self.execute(cmd)
        except Exception as err:
            raise(err)
        finally:
            # clean symlink
            os.unlink(fname)
