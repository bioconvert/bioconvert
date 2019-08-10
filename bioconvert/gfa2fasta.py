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

"""Convert :term:`GFA` to :term:`FASTA` format"""

import colorlog

from bioconvert.core.decorators import requires, requires_nothing, compressor
from bioconvert import ConvBase

logger = colorlog.getLogger(__name__)


__all__ = ["GFA2FASTA"]


class GFA2FASTA(ConvBase):
    """Convert sorted :term:`GFA` file into :term:`FASTA` file 

    Available methods are based on awk or python (default)

    .. plot::

         from bioconvert.gfa2fasta import GFA2FASTA
         from bioconvert import bioconvert_data
         from easydev import TempFile

         with TempFile(suffix=".fasta") as fh:
             infile = bioconvert_data("test_gfa2fasta_v1.gfa")
             convert = GFA2FASTA(infile, fh.name)
             convert.boxplot_benchmark()

    :reference: https://github.com/GFA-spec/GFA-spec/blob/master/GFA-spec.md

    .. seealso:: :mod:`bioconvert.simulator.gfa`

    """
    _default_method = "python"

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input BAM file. **It must be sorted**.
        :param str outfile: The path to the output file
        """
        super().__init__(infile, outfile)

    @requires("awk")
    @compressor
    def _method_awk(self, *args, **kwargs):
        """

        :return: the standard output
        :rtype: :class:`io.StringIO` object.

        .. note:: this method fold the sequence to 80 characters
        """
        # Note1: since we use .format, we need to escape the { and } characters
        # Note2: the \n need to be escaped for Popen to work
        cmd = """awk '/^S/{{print ">"$2"\\n"$3}}' {} | fold > {}""".format(self.infile, self.outfile)
        self.execute(cmd)

    @requires_nothing
    @compressor
    def _method_python(self, *args, **kwargs):
        with open(self.infile, "r") as fin:
            with open(self.outfile, "w") as fout:
                for i, line in enumerate(fin.readlines()):
                    if line.startswith("S"):
                        args = line.split()
                        if len(args) == 3:
                            fout.write(">{}\n{}\n".format(args[1], args[2]))
                        elif len(args) == 4:
                            fout.write(">{}\n{}\n".format(args[1]+" " + args[3], args[2]))
                        else:
                            raise ValueError("Illformed line on line {}. Expected 3 or 4 values".format(i))



