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

"""Convert :term:`MAF` file to :term:`SAM` format"""
from bioconvert import ConvBase

import colorlog
logger = colorlog.getLogger(__name__)

__all__ = ["MAF2SAM"]


class MAF2SAM(ConvBase):
    """This is the Multiple alignment format or MIRA assembly format

    This is not Mutation Annotation Format (somatic)


    pbsim creates this kind of data

    Some references:

    - https://github.com/peterjc/maf2sam/
    - https://github.com/arq5x/nanopore-scripts/master/maf-convert.py
    - http://bioperl.org/formats/alignment_formats/MAF_multiple_alignment_format.html

    Those two codes were in Py2 at the time of this implementation. We re-used
    some of the information from maf-convert but the code in
    bioconvert.io.maf can be considered original. 
    """

    def __init__(self, infile, outfile):
        super().__init__(infile, outfile)

    def _method_python(self, *args, **kwargs):
        from bioconvert.io import maf
        conv = maf.MAF(self.infile, self.outfile)
        conv.to_sam()

