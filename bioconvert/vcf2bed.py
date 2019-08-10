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

"""Convert :term:`VCF`  to :term:`BED3` file"""
from bioconvert import ConvBase
import colorlog

from bioconvert.core.decorators import requires

logger = colorlog.getLogger(__name__)


__all__ = ["VCF2BED"]


class VCF2BED(ConvBase):
    """
    Convert VCF file to BED3 file by extracting positions.

    The awk method implemented here below reports an interval
    of 1 for SNP, the length of the insertion or the length of
    the deleted part in case of deletion.
    """
    _default_method = "awk"

    @requires("awk")
    def _method_awk(self, *args, **kwargs):
        """
        do the conversion :term:`VCF` -> :term:`BED` using awk

        :return: the standard output
        :rtype: :class:`io.StringIO` object.
        """
        awkcmd = """awk '{{if(length($4) > length($5)) print $1,($2-1),($2+length($4)-1); else print $1,($2-1),($2+length($5)-1)}}' OFS='\t'"""
        cmd = """awk '! /\#/' {} | {} > {}""".format(self.infile, awkcmd, self.outfile)
        self.execute(cmd)
