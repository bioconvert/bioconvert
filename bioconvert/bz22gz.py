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
"""Convert :term:`BZ2` to :term:`GZ` format"""
import bz2
import gzip

import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires, requires_nothing

logger = colorlog.getLogger(__name__)


__all__ = ["BZ22GZ"]


class BZ22GZ(ConvBase):
    """Convert :term:`BZ2` file to :term:`GZ` file

    Methods based on bunzip2 or zlib/bz2 Python libraries.

    """

    #: Default value
    _default_method = "bz2_gz"
    _threading = True

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input BZ2 file
        :param str outfile: output GZ filename

        """
        super(BZ22GZ, self).__init__(infile, outfile, *args, **kargs)

    @requires("bunzip2")
    def _method_bz2_gz(self, *args, **kwargs):
        """Method that uses bunzip2 gzip.

        `bunzip2 documentation <https://docs.oracle.com/cd/E86824_01/html/E54763/bunzip2-1.html>`_
        `gzip documentation <https://www.gnu.org/software/gzip/manual/gzip.html>`_"""
        # conversion
        cmd = "bunzip2 -c {input} | gzip > {output}".format(input=self.infile, output=self.outfile)
        self.execute(cmd)

    @requires_nothing
    def _method_python(self, *args, **kargs):
        """Internal method"""
        with bz2.open(self.infile, "rb") as f, gzip.open(self.outfile, "wb") as g:
            g.write(f.read())
