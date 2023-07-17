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
"""Convert :term:`FASTQ` to :term:`FASTA` format"""
from bioconvert import ConvBase, bioconvert_script

from bioconvert.core.decorators import requires, requires_nothing


class FAST52POD5(ConvBase):
    """Convert :term:`FAST5` to :term:`POD5`


    """

    _default_method = "pod5"
    _loss = False
    _threading = True  # not used but required for the main benchmark

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input FAST5 file.
        :param str outfile: The path to the output POD5 file.
        """
        super(FAST52POD5, self).__init__(infile, outfile)

    @requires("pod5")
    def _method_pod5(self, *args, **kwargs):
        """
        """
        if kwargs.get('force', False) is True:
            cmd = f"pod5 convert fast5 {self.infile} --output {self.outfile} --force-overwrite "
        else:
            cmd = f"pod5 convert fast5 {self.infile} --output {self.outfile}"
        self.execute(cmd)

