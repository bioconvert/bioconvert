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
"""Convert :term:`BIGBED` format to :term:`BED` format """
import os

import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


__all__ = ["BIGBED2BED"]


class BIGBED2BED(ConvBase):
    """Converts a sequence alignment in :term:`BIGBED` format to :term:`BED4` format

    Methods available are based on pybigwig [DEEPTOOLS]_.
    """
    _default_method = 'pybigwig'

    def __init__(self, infile, outfile):#=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`BIGBED` file.
        :param str outfile: (optional) output :term:`BED4` file
        """
        super(BIGBED2BED, self).__init__(infile, outfile)

    @requires(python_library="pyBigWig")
    def _method_pybigwig(self, *args, **kwargs):
        """
        """
        import pyBigWig
        bw = pyBigWig.open(self.infile)
        assert bw.isBigBed() is True, "Not a valid bigBed file"

        with open(self.outfile, "w") as fout:

            for chrom in  sorted(bw.chroms().keys()):
                L = bw.chroms()[chrom]
                for tup in bw.entries(chrom, 0, L):
                    s, e, val = tup
                    # the bigbed may have many columns. We consider only bed4
                    # for now
                    # FIXME
                    try:
                        val = val.split("\t")[0]
                    except:
                        pass
                    fout.write("{}\t{}\t{}\t{}\n".format(chrom, s, e, val))
