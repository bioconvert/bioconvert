# -*- coding: utf-8 -*-

###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018  Institut Pasteur, Paris and CNRS.                     #
# See the COPYRIGHT file for details                                      #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
###########################################################################

"""Convert :term:`SAM` file to :term:`CRAM` format"""
import os

from bioconvert import ConvBase
from bioconvert.core.base import ConvArg

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)

__all__ = ["SAM2CRAM"]


class SAM2CRAM(ConvBase):
    """Convert :term:`SAM` file to :term:`CRAM` file

    The conversion requires the reference corresponding to the input file
    It can be provided as an argument with the standalone (*-\\-reference*). 
    Otherwise, users are asked to provide it.

    Methods available are based on samtools [SAMTOOLS]_.
    """
    _default_method = "samtools"
    _threading = True

    def __init__(self, infile, outfile, reference=None, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input SAM file
        :param str outfile: output CRAM filename

        """
        super(SAM2CRAM, self).__init__(infile, outfile, *args, **kargs)

    def _get_reference(self):
        # In case the --reference is not used
        msg = "Please enter the reference corresponding "
        msg += "to the input BAM file:"
        reference = input(msg)
        if os.path.exists(reference) is False:
            raise IOError("Reference required")
        else:
            logger.debug("Reference exists ({}).".format(reference))
        return reference

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        # -C means output is CRAM

        reference = kwargs.get("reference", None)
        if reference is None:
            reference = self._get_reference()

        cmd = "samtools view -@ {} -C {} -T {} > {}".format(
            self.threads, self.infile, reference, self.outfile)
        try:
            self.execute(cmd)
        except:
            logger.debug("FIXME. The ouput message from samtools is on stderr...")

    @classmethod
    def get_additional_arguments(cls):
        yield ConvArg(
            names="--reference",
            nargs=1,
            default=None,
            #type=ConvArg.file,
            help="reference used",
        )
