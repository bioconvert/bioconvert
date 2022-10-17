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
"""Convert :term:`FASTA` to :term:`PHYLIP` format"""
import colorlog
from Bio import SeqIO

from bioconvert import ConvBase
from bioconvert.core.decorators import compressor, requires

_log = colorlog.getLogger(__name__)


__all__ = ["FASTA2PHYLIP"]


class FASTA2PHYLIP(ConvBase):
    """Converts a sequence alignment in :term:`FASTA` format to :term:`PHYLIP` format

    Conversion is based on Bio Python modules

    Methods available are based on squizz [SQUIZZ]_ or biopython [BIOPYTHON]_ or
    goalign [GOALIGN]_. Squizz is the default
    (https://github.com/bioconvert/bioconvert/issues/149). Phylip created is a
    strict phylip that is with 10 characters on the first column.

    """

    #: Default value
    _default_method = "biopython"

    def __init__(self, infile, outfile=None, alphabet=None, *args, **kwargs):
        """.. rubric:: constructor

        :param str infile: input :term:`FASTA` file.
        :param str outfile: (optional) output :term:`PHYLIP` file
        """
        super(FASTA2PHYLIP, self).__init__(infile, outfile)

    @requires(python_library="biopython")
    @compressor
    def _method_biopython(self, *args, **kwargs):
        """For this method we use the biopython package Bio.SeqIO.

        `Bio.SeqIO Documentation <https://biopython.org/docs/1.76/api/Bio.SeqIO.html>`_"""
        sequences = list(SeqIO.parse(self.infile, "fasta"))
        count = SeqIO.write(sequences, self.outfile, "phylip")
        _log.debug("Converted %d records to phylip" % count)

    @requires("squizz")
    @compressor
    def _method_squizz(self, *args, **kwargs):
        """Convert fasta file in Phylip interleaved format using squizz tool.
        The fasta file must be an alignement file, this means that all sequences must
        have the same length (with the gap) otherwise an error will be raised."""
        cmd = "squizz -c PHYLIPI {infile} > {outfile}".format(infile=self.infile, outfile=self.outfile)
        self.execute(cmd)

    @requires("go")
    @compressor
    def _method_goalign(self, *args, **kwargs):
        """Convert fasta file in Phylip interleaved format using goalign tool.

        `goalign documentation <https://github.com/fredericlemoine/goalign>`_

        The fasta file must be an alignemnt file, this means that  all sequences
        must have the same length (with the gap) otherwise an error will be raised
        """
        self.install_tool("goalign")
        cmd = "goalign reformat phylip -i {infile} -o {outfile}".format(infile=self.infile, outfile=self.outfile)
        self.execute(cmd)
