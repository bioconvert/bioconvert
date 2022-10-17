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
"""Convert :term:`FASTQ` to :term:`QUAL` format"""
from bioconvert import ConvBase, bioconvert_script, logger
from bioconvert.core.base import ConvArg
from bioconvert.core.decorators import compressor, in_gz, out_compressor, requires, requires_nothing

logger.__name__ = "fastq2qual"


class FASTQ2QUAL(ConvBase):
    """Convert :term:`FASTQ` to :term:`QUAL`"""

    # use readfq for now because pure python are fast enough
    # for production, could use seqtk which seems the fastest method though
    # Make sure that the default handles also the compresssion
    #: Default value
    _default_method = "readfq"

    # (https://raw.githubusercontent.com/lh3/readfq/master/readfq.py)
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
                if l[0] in "@+":
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield header, seq, "".join(seqs)  # yield a fastq record
                    break

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FASTA file.
        :param str outfile: The path to the output file.
        """
        super(FASTQ2QUAL, self).__init__(infile, outfile)

    @requires_nothing
    @compressor
    def _method_readfq(self, *args, **kwargs):
        """This method is inspired by Readfq coded by Heng Li.

        `original Readfq method <https://github.com/lh3/readfq>`_"""
        with open(self.outfile, "w") as outfile, open(self.infile, "r") as fastq:
            for (name, seq, qual) in FASTQ2QUAL._readfq(fastq):
                outfile.write(">{}\n{}\n".format(name, qual))
