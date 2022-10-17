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
"""Convert :term:`GENBANK` to :term:`GFF3` format"""

from bioconvert import ConvBase

__all__ = ["GENBANK2GFF3"]


class GENBANK2GFF3(ConvBase):
    """Convert :term:`GENBANK` file to :term:`GFF3` file

    Method based on biocode.

    """

    #: Default value
    _default_method = "biocode"

    def __init__(self, infile, outfile, *args, **kargs):
        """.. rubric:: constructor

        :param str infile: input GENBANK file
        :param str outfile: output GFF3 filename

        """
        super(GENBANK2GFF3, self).__init__(infile, outfile)

    def _method_biocode(self, *args, **kwargs):
        """Uses scripts from biocode copied and modified in bioconvert.utils.biocode

        Please see `Main entry  <https://github.com/jorvis/biocode/>`_
        """
        from bioconvert.utils.biocode.convert_genbank_to_gff3 import gbk2gff3

        gbk2gff3(self.infile, self.outfile, fasta=False)
