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
"""Convert :term:`BCF` file to :term:`VCF` format"""
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


class BCF2VCF(ConvBase):
    """Convert :term:`BCF` file to :term:`VCF` file

    Methods available are based on bcftools [BCFTOOLS]_.

    """

    _default_method = "bcftools"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input BCF file
        :param str outfile: output VCF file

        """
        super(BCF2VCF, self).__init__(infile, outfile, *args, **kargs)

    @requires("bcftools")
    def _method_bcftools(self, *args, **kwargs):
        """Here we use the bcftools tool from samtools.

        `bcftools documentation <https://samtools.github.io/bcftools/bcftools.html>`_"""

        # -O, --output-type b|u|z|v Output compressed BCF (b), uncompressed BCF
        # (u), compressed VCF (z), uncompressed VCF (v). Use the -Ou option when
        # piping between bcftools subcommands to speed up performance
        cmd = "bcftools view {} -O v -o {}".format(self.infile, self.outfile)
        self.execute(cmd)
