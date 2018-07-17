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

from bioconvert import ConvBase
from bioconvert.core.decorators import requires_nothing
from bioconvert.readers.fasta import Fasta

__all__ = ["FASTA2BED"]


class FASTA2BED(ConvBase):
    """Convert :term:`FASTA` file to :term:`BED` file
    This convertion is not straight forward, so, please refer to method documentations
    for precise convertion details
    """
    _default_method = "scaffold"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input FASTA file
        :param str outfile: output BED filename

        """
        super(FASTA2BED, self).__init__(infile, outfile, *args, **kargs)

    @requires_nothing
    def _method_scaffold(self, *args, **kwargs):
        """ This method will convert a fasta into a bed file assuming that the fasta represent a scaffold
        Each sequence in the file will be considered as following the previous one in the scaffold.
        For this pseudo-annotation, the contigs in the scaffold will be considered as neighbors without gaps.
        """
        reader = Fasta(self.infile)

        previous_idx = 0
        with open(self.outfile, "w") as writer:
            for sequence in reader.read():
                writer.write("{}\t{}\t{}\n".format(sequence["id"], previous_idx, previous_idx+len(sequence["value"])-1))
                previous_idx += len(sequence["value"])
