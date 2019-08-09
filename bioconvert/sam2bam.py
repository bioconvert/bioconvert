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

"""Convert :term:`SAM` file to :term:`BAM` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

import colorlog

logger = colorlog.getLogger(__name__)

__all__ = ["SAM2BAM"]


class SAM2BAM(ConvBase):
    """Convert :term:`SAM` file to :term:`BAM` file"""

    _threading = True

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:

        """
        super(SAM2BAM, self).__init__(infile, outfile, *args, **kargs)

    @requires("samtools")
    def _method_samtools(self, *args, **kwargs):
        """ Do the conversion :term:`SAM` -> :term:`BAM` using samtools"""
        # -S means ignored (input format is auto-detected)
        # -b means output to BAM format
        # -h means include header in SAM output
        cmd = "samtools view -Sbh -@ {} {} > {}".format(self.threads,
                                                        self.infile,
                                                        self.outfile)
        self.execute(cmd)
