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
"""Convert :term:`GZ` to :term:`DSRC` format"""

import colorlog

from bioconvert import ConvBase
from bioconvert import requires

logger = colorlog.getLogger(__name__)


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
    _default_method = "pigzdsrc"
    _threading = True

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GZ filename
        :param str outfile: output DSRC filename

        """
        super(GZ2DSRC, self).__init__(infile, outfile, *args, **kargs)

    @requires(external_binaries=["pigz", "dsrc"])
    def _method_pigzdsrc(self, *args, **kwargs):
        """
        do the conversion gz -> :term:`DSRC`

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        cmd = "pigz -d -c -p {threads} {input} | dsrc c -s -t{threads} {output}"
        self.execute(cmd.format(
            threads=self.threads,
            input=self.infile,
            output=self.outfile))


