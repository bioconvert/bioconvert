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

"""Convert :term:`NEXUS` to :term:`FASTA` format"""

import colorlog

from bioconvert import ConvBase
from bioconvert.core.decorators import requires
from bioconvert.core.decorators import compressor

_log = colorlog.getLogger(__name__)


__all__ = ['NEXUS2FASTA']


class NEXUS2FASTA(ConvBase):
    """
    Converts a sequence alignment from :term:`NEXUS` format to :term:`FASTA` format.
    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`NEXUS` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires("go")
    @compressor
    def _method_goalign(self, *args, **kwargs):
        """
        Convert :term:`NEXUS` interleaved file in :term:`FASTA` format using goalign tool.
        https://github.com/fredericlemoine/goalign

        .. warning::
            the sequential format is not supported

        """
        self.install_tool('goalign')
        cmd = 'goalign reformat fasta x -i {infile} -o {outfile} -x'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """
        Convert :term:`NEXUS` interleaved or sequential file in :term:`FASTA` format using biopython.
        The FASTA output file will be an aligned FASTA file

For instance:

We have a Nexus input file that look like ::

    #NEXUS
    [TITLE: Test file]

    begin data;
    dimensions ntax=3 nchar=123;
    format interleave datatype=DNA missing=N gap=-;

    matrix
    read3                -AT--------CCCGCTCGATGGGCCTCATTGCGTCCACTAGTTGATCTT
    read2                -----------------------GGAAGCCCACGCCACGGTCTTGATACG
    read4                ---------------------AGGGATGAACGATGCTCGCAGTTGATGCT

    read3                CTGGAGTAT---T----TAGGAAAGCAAGTAAACTCCTTGTACAAATAAA
    read2                AATTTTTCTAATGGCTATCCCTACATAACCTAACCGGGCATGTAATGTGT
    read4                CAGAAGTGCCATTGCGGTAGAAACAAATGTTCCCAGATTGTTGACTGATA

    read3                GATCTTA-----GATGGGCAT--
    read2                CACCGTTGTTTCGACGTAAAGAG
    read4                AGTAGGACCTCAGTCGTGACT--
    ;

    end;
    begin assumptions;
    options deftype=unord;
    end;

the output file will look like ::

    >read3
    -AT--------CCCGCTCGATGGGCCTCATTGCGTCCACTAGTTGATCTTCTGGAGTAT-
    --T----TAGGAAAGCAAGTAAACTCCTTGTACAAATAAAGATCTTA-----GATGGGCA
    T--
    >read2
    -----------------------GGAAGCCCACGCCACGGTCTTGATACGAATTTTTCTA
    ATGGCTATCCCTACATAACCTAACCGGGCATGTAATGTGTCACCGTTGTTTCGACGTAAA
    GAG
    >read4
    ---------------------AGGGATGAACGATGCTCGCAGTTGATGCTCAGAAGTGCC
    ATTGCGGTAGAAACAAATGTTCCCAGATTGTTGACTGATAAGTAGGACCTCAGTCGTGAC
    T--

and not ::

    >read3
    ATCCCGCTCGATGGGCCTCATTGCGTCCACTAGTTGATCTTCTGGAGTATTTAGGAAAGC
    AAGTAAACTCCTTGTACAAATAAAGATCTTAGATGGGCAT
    >read2
    GGAAGCCCACGCCACGGTCTTGATACGAATTTTTCTAATGGCTATCCCTACATAACCTAA
    CCGGGCATGTAATGTGTCACCGTTGTTTCGACGTAAAGAG
    >read4
    AGGGATGAACGATGCTCGCAGTTGATGCTCAGAAGTGCCATTGCGGTAGAAACAAATGTT
    CCCAGATTGTTGACTGATAAGTAGGACCTCAGTCGTGACT
"""
        from Bio import AlignIO
        with open(self.outfile, "w") as output_handle:
            alignments = list(AlignIO.parse(self.infile, "nexus", alphabet=self.alphabet))
            AlignIO.write(alignments, output_handle, "fasta")

    @requires("squizz")
    @compressor
    def _method_squizz(self, *args, **kwargs):
        """
        Convert :term:`NEXUS` sequential or interleave file in :term:`FASTA` format using squizz tool.

        command used::

            squizz -c FASTA infile > outfile
        """
        cmd = 'squizz -c FASTA {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
