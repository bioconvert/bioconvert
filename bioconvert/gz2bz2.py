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
"""Convert :term:`GZ` file to :term:`BZ2` format"""
import bz2
import gzip

from bioconvert import ConvBase
from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires, requires_nothing

__all__ = ["GZ2BZ2"]


class GZ2BZ2(ConvBase):
    """Convert :term:`GZ` file to :term:`BZ2` file

    Unzip input file using pigz or gunzip and compress using pbzip2. Default
    is pigz/pbzip2.

    """

    _threading = True

    #: Default value
    _default_method = "pigz_pbzip2"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GZ file
        :param str outfile: output BZ2 filename

        """
        super(GZ2BZ2, self).__init__(infile, outfile, *args, **kargs)

    @requires(
        external_binaries=[
            "pigz",
            "pbzip2",
        ]
    )
    def _method_pigz_pbzip2(self, *args, **kwargs):
        """Method that uses pigz pbzip2.

        `pigz documentation <https://linux.die.net/man/1/pigz>`_
        `pbzip2 documentation <https://linux.die.net/man/1/pbzip2>`_"""
        # conversion
        cmd = "pigz -d -c -p {threads} {input} | pbzip2 -p{threads} > {output}"
        self.execute(cmd.format(threads=self.threads, input=self.infile, output=self.outfile))

    @requires(
        external_binaries=[
            "gunzip",
            "bzip2",
        ]
    )
    def _method_gunzip_bzip2(self, *args, **kwargs):
        """Single theaded conversion. Method that uses gunzip bzip2.

        `gunzip documentation <https://linux.die.net/man/1/gunzip>`_
        `bzip2 documentation <http://www.delafond.org/traducmanfr/man/man1/bzip2.1.html>`_"""
        cmd = "gunzip --to-stdout {input} | bzip2 > {output}"
        self.execute(cmd.format(input=self.infile, output=self.outfile))

    @requires_nothing
    def _method_python(self, *args, **kwargs):
        """Internal method"""
        with gzip.open(self.infile, "rb") as f, bz2.open(self.outfile, "wb") as g:
            g.write(f.read())
