###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
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
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
""" GFF3 to GTF conversion """

from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires, requires_nothing
from bioconvert import ConvBase

__all__ = ["GFF32GTF"]


class GFF32GTF(ConvBase):
    """Convert :term:`GFF3` file to :term:`GTF` file"""

    _default_method = "gffread"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GFF3 filename
        :param str outfile: output GTF filename

        """
        super(GFF32GTF, self).__init__(infile, outfile, *args, **kargs)

    @requires("gffread")
    def _method_gffread(self, *args, **kwargs):
        """some description"""
        cmd = f"gffread -T {self.infile} -o {self.outfile}"
        # use self.infile, self.outfile
        self.execute(cmd)
