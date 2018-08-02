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
import os
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import requires

_log = colorlog.getLogger(__name__)


class CLUSTAL2FASTA(ConvBase):
    """
    Converts a sequence alignment from :term:`CLUSTAL` format to :term:`FASTA` format. ::

        converter = CLUSTAL2FASTA(infile, outfile)
        converter(method='biopython')

    default method = biopython
    available methods = biopython, squizz, goalign
    """
    _default_method = 'biopython'

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`CLUSTAL` file.
        :param str outfile: (optional) output :term:`FASTA` file
        """
        super().__init__(infile, outfile)
        self.alphabet = alphabet

    @requires(python_library="biopython")
    def _method_biopython(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` interleaved file in :term:`PHYLIP` format using biopython.

        :param threads: not used.
        """
        sequences = list(SeqIO.parse(self.infile, "clustal", alphabet=self.alphabet))
        count = SeqIO.write(sequences, self.outfile, "fasta")
        _log.info("Converted %d records to fasta" % count)

    @requires("squizz")
    def _method_squizz(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` file in :term:`FASTA` format using squizz tool.

        :param threads: not used.
        """
        cmd = 'squizz -c FASTA {infile} > {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)

    @requires("conda")
    def _method_goalign(self, threads=None, *args, **kwargs):
        """
        Convert :term:`CLUSTAL` file in  :term:`FASTA` format using goalign tool.
        https://github.com/fredericlemoine/goalign

        :param threads: not used.
        """
        self.install_tool('goalign')
        cmd = 'goalign reformat fasta --clustal -i {infile} -o {outfile}'.format(
            infile=self.infile,
            outfile=self.outfile)
        self.execute(cmd)
