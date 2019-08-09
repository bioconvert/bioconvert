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
"""Convert :term:`BAM` format to :term:`JSON` format"""
from bioconvert import ConvBase

import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


class BAM2JSON(ConvBase):
    """Convert :term:`BAM` format to :term:`JSON` file

    Methods available are based on bamtools [BAMTOOLS]_.

    """

    _default_method = "bamtools"

    def __init__(self, infile, outfile):
        """.. rubric:: constructor

        :param str infile:
        :param str outfile:
        """
        super(BAM2JSON, self).__init__(infile, outfile)

    @requires("bamtools")
    def _method_bamtools(self, *args, **kwargs):
        """Do the conversion :term:`BAM` -> :term:`JSON` using bamtools """

        cmd = "bamtools convert -format json -in {0} -out {1}".format(
            self.infile, self.outfile
        )
        self.execute(cmd)
