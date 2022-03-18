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
"""Convert :term:`VCF` format to :term:`PLINK` formats"""

import os
import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.utils import generate_outfile_name

_log = colorlog.getLogger(__name__)


class VCF2PLINK(ConvBase):
    """Converts a genotype dataset vcf :term:`VCF` format to
    ped+map in :term:`PLINK` format

    Conversion is based on plink executable

    """

    #: Default value
    _default_method = "plink"

    def __init__(self, infile, outfile=None, *args, **kwargs):
        """For this method, we use the plink tool.
        
        `plink documentation <http://hpc.ilri.cgiar.org/beca/training/data_mgt_2017/BackgroundMaterial/PlinkTutorial.pdf>`_
        
        .. rubric:: constructor

        :param str infile: input :term:`PLINK` file.
        :param str outfile: (optional) output :term:`BPLINK` file
        """
        if not outfile:
            outfile = os.path.splitext(infile)[0]
        super(VCF2PLINK, self).__init__(infile, outfile)

    @requires("plink")
    def _method_plink(self, *args, **kwargs):
        """
        Convert using plink executable.
        """
        cmd = "plink --vcf {infile} --recode --out {outfile}".format(
            infile=self.infile, outfile=self.outfile
        )
        self.execute(cmd)
