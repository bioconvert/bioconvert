
# -*- coding: utf-8 -*-
###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Authors: see CONTRIBUTORS.rst                                           #
# Copyright © 2018-2019  Institut Pasteur, Paris and CNRS.                #
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
"""Convert :term:`WIG` to :term:`BED` format"""

from bioconvert import ConvBase
from bioconvert import requires, compressor

__all__ = ["WIG2BED"]


class WIG2BED(ConvBase):
    """Convert :term:`WIG` file to :term:`BED` file

    Method availabe are based on wig2bed tool.

    """

    _default_method = "bedops"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input WIG file
        :param str outfile: output BED filename

        """
        super(WIG2BED, self).__init__(infile, outfile, *args, **kargs)

    @requires("wig2bed")
    @compressor
    def _method_default(self, *args, **kwargs):
        """some description"""
        cmd = "wig2bed < {} > {}".format(self.infile, self.outfile)
        self.execute(cmd)

