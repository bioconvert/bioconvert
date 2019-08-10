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
"""Convert :term:`SCF` file to :term:`FASTQ` file"""
import sys
import struct

from collections import defaultdict
from bioconvert import ConvBase
from bioconvert.io import scf
import colorlog

from bioconvert.core.decorators import requires_nothing, compressor

_log = colorlog.getLogger(__name__)

__all__ = ["SCF2FASTQ"]


class SCF2FASTQ(ConvBase):
    """
    Converts a binary :term:`SCF` file to :term:`FastQ` file

    :param str infile: input SCF file
    :param str outfile: output name file
    """

    @requires_nothing
    @compressor
    def _method_python(self, *args, **kwargs):
        sequence, qualities, comments = scf.read_scf(self.infile)

        # Wrinting output file
        with open(self.outfile, "w") as output_file:
            output_file.write("@" + comments.replace("\n", "-").replace(" ", "_") + "\n")
            output_file.write(sequence + "\n")
            output_file.write("+" + comments.replace("\n", "-").replace(" ", "_") + "\n")
            for i in qualities:
                if i > 92:
                    output_file.write(chr(126))
                else:
                    output_file.write(chr(i+34))
            output_file.write("\n")

        """
        print(sequence)
        print("INFORMATIONS")
        print("magic_number = " + magic_number.decode("utf-8"))
        print("samples = " + str(samples))
        print("samples_offset = " + str(samples_offset))
        print("bases = " + str(bases))
        print("bases_left_clip = " + str(bases_left_clip))
        print("bases_right_clip = " + str(bases_right_clip))
        print("bases_offset = " + str(bases_offset))
        print("comments_size = " + str(comments_size))
        print("comments_offset = " + str(comments_offset))
        print("version = " + version.decode("utf-8"))
        print("sample_size = " + str(sample_size))
        print("code_set = " + str(code_set))
        print("private_size = " + str(private_size))
        print("private_offset = " + str(private_offset))
        #print("spare = " + str(spare))
        print()
        print("COMMENTS")
        print(comments)
        print()
        print("SEQUENCE")
        print(sequence)
        print()
        print("QUALITIES")
        print(qualities)
        """




"""
http://staden.sourceforge.net/manual/formats_unix_2.html
http://doc.bioperl.org/bioperl-live/Bio/SeqIO/scf.html#POD6
https://docs.python.org/2/library/struct.html
http://www.perlmonks.org/?node_id=224666

SCF file organisation (more or less)

Length in bytes                        Data
---------------------------------------------------------------------------
128                                    header
Number of samples * sample size        Samples for A trace
Number of samples * sample size        Samples for C trace
Number of samples * sample size        Samples for G trace
Number of samples * sample size        Samples for T trace
Number of bases * 4                    Offset into peak index for each base
Number of bases                        Accuracy estimate bases being 'A'
Number of bases                        Accuracy estimate bases being 'C'
Number of bases                        Accuracy estimate bases being 'G'
Number of bases                        Accuracy estimate bases being 'T'
Number of bases                        The called bases
Number of bases * 3                    Reserved for future use
Comments size                          Comments
Private data size                      Private data
"""

