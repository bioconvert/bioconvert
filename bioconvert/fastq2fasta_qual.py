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
"""Convert :term:`FASTQ` to :term:`FASTA` and :term:`QUAL` formats"""
from bioconvert import ConvBase, bioconvert_script
from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import compressor, in_gz
from bioconvert.core.decorators import requires, requires_nothing

from mappy import fastx_read
import mmap


__all__ = ["FASTQ2FASTA_QUAL"]


class FASTQ2FASTA_QUAL(ConvBase):
    """Convert :term:`FASTQ` to :term:`FASTA` and :term:`QUAL` formats

    Method based on Python.

    """
    _default_method = "python"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input fastQ file
        :param str outfile: output fastA filename
        :param str outfile2: output QUAL file

        """
        super(FASTQ2FASTA_QUAL, self).__init__(infile, outfile)

        # Here we need to reset the infile and outfile since we have
        # a 1-to-2 conversion.
        self.infile = infile
        self.outfile = outfile[0]
        self.outfile2 = outfile[1]

    @staticmethod
    def _readfq(fp):  # this is a generator function
        last = None  # this is a buffer keeping the last unprocessed line
        while True:  # mimic closure; is it a bad idea?
            if not last:  # the first record or a record following a fastq
                for l in fp:  # search for the start of the next record
                    if l[0] == "@":  # fastq header line
                        last = l[:-1]  # save this line
                        break
            if not last:
                break
            header, seqs, last = last[1:], [], None
            for l in fp:  # read the sequence
                if l[0] in '@+':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield header, seq, ''.join(seqs)  # yield a fastq record
                    break

    @requires_nothing
    @compressor
    def _method_python(self, *args, **kwargs):
        with open(self.outfile, "w") as fasta, open(self.outfile2, "w") as quality, open(self.infile, "r") as fastq:
            for (name, seq, qual) in FASTQ2FASTA_QUAL._readfq(fastq):
                fasta.write(">{}\n{}\n".format(name, seq))
                quality.write(">{}\n{}\n".format(name, qual))

    @staticmethod
    def get_IO_arguments():
        yield ConvArg(
            names="input_file",
            default=None,
            type=ConvArg.file,
            help="The path to the file to convert.",
        )
        yield ConvArg(
            names="output_file",
            nargs= 2,
            default=None,
            type=ConvArg.file,
            output_argument=True,
            help="The path where the result will be stored.",
        )
