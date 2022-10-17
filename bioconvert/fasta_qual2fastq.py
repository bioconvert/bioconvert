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
"""Convert :term:`FASTA` format to :term:`FASTQ` format"""
import sys

import colorlog

from bioconvert import ConvBase
from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import requires
from bioconvert.core.extensions import extensions

_log = colorlog.getLogger(__name__)


class FASTA_QUAL2FASTQ(ConvBase):
    """Convert FASTA and QUAL back into a FASTQ file

    Method based on pysam [PYSAM]_.

    """

    #: Default value
    _default_method = "pysam"

    def __init__(self, infile, outfile):
        """
        :param list infile: The path to the input FASTA file, the path to the input QUAL file
        :param str outfile: The path to the output FASTQ file
        """
        super(FASTA_QUAL2FASTQ, self).__init__(infile, outfile)

    @requires(python_library="pysam")
    def _method_pysam(self, *args, **kwargs):
        """This method uses the FastxFile function of the Pysam python module.

        `FastxFile documentation <https://pysam.readthedocs.io/en/latest/api.html#pysam.FastxFile.close>`_"""
        from pysam import FastxFile

        if self.infile[1] is None:
            _log.error("No quality file provided. Please add a quality file path ")
            sys.exit(1)

        else:  # length must be equal and identifiers sorted similarly
            with open(self.outfile, "w") as fastq_out:
                for seq, qual in zip(FastxFile(self.infile[0]), FastxFile(self.infile[1])):
                    assert seq.name == qual.name
                    if seq.comment:
                        fastq_out.write(
                            "@{0} {1}\n{2}\n+\n{3}\n".format(seq.name, seq.comment, seq.sequence, qual.sequence)
                        )
                    else:
                        fastq_out.write("@{0}\n{1}\n+\n{2}\n".format(seq.name, seq.sequence, qual.sequence))

    @staticmethod
    def get_IO_arguments():
        yield ConvArg(
            names="input_file",
            nargs=2,
            default=None,
            type=ConvArg.file,
            help="The path to the file to convert.",
        )
        yield ConvArg(
            names="output_file",
            default=None,
            type=ConvArg.file,
            output_argument=True,
            help="The path where the result will be stored.",
        )
