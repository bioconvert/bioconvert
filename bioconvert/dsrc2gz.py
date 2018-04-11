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
""" Convert a compressed fastq.gz file to :term:`DSRC` compression format """
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


__all__ = ["DSRC2GZ"]


class DSRC2GZ(ConvBase):
    """Convert compressed fastq.dsrc file into fastq.gz compressed file

    .. plot::

         from bioconvert.dsrc2gz import DSRC2GZ
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".gz") as fh:
             infile = bioconvert_data("test_SP1.fq.dsrc")
             convert = DSRC2GZ(infile, fh.name)
             convert.boxplot_benchmark()

    """
    _default_method = "dsrcpigz"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input DSRC filename
        :param str outfile: output GZ filename

        """
        super(DSRC2GZ, self).__init__(infile, outfile, *args, **kargs)

    @requires("dsrc")
    def _method_dsrcpigz(self, *args, **kwargs):
        """
        do the conversion dsrc -> :term:'GZ`
        """
        cmd = "dsrc d -s -t {threads} {input} | pigz -c -p {threads} > {output}"
        self.execute(cmd.format(
            threads=self.threads,
            input=self.infile,
            output=self.outfile))


