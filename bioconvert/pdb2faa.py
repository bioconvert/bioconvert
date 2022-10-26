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
"""Convert :term:`PDB` to :term:`FAA` format"""
import sys, re
from bioconvert import ConvBase

from bioconvert.core.decorators import requires_nothing


class PDB2FAA(ConvBase):
    """Convert :term:`PDB` to :term:`FAA`

    Convert PDB file into amino acid file
    """

    _default_method = "bioconvert"

    def __init__(self, infile, outfile):
        """
        :param str infile: The path to the input PDB file.
        :param str outfile: The path to the output file.
        """
        super(PDB2FAA, self).__init__(infile, outfile)

    @requires_nothing
    def _method_bioconvert(self, *args, **kwargs):

        long2short = {
            'ALA':'A', 'VAL':'V', 'PHE':'F', 'PRO':'P', 'MET':'M',
           'ILE':'I', 'LEU':'L', 'ASP':'D', 'GLU':'E', 'LYS':'K',
           'ARG':'R', 'SER':'S', 'THR':'T', 'TYR':'Y', 'HIS':'H',
           'CYS':'C', 'ASN':'N', 'GLN':'Q', 'TRP':'W', 'GLY':'G',
           'MSE':'M',
        }

        pattern=re.compile("^ATOM\s{2,6}\d{1,5}\s{2}CA\s[\sA]([A-Z]{3})\s([\s\w])|^HETATM\s{0,4}\d{1,5}\s{2}CA\s[\sA](MSE)\s([\s\w])")

        chain_dict=dict()
        chain_list=[]

        with open(self.infile, 'r') as fin:
            for line in fin.readlines():
                if line.startswith("ENDMDL"):
                    break
                match_list = pattern.findall(line)
                if match_list:
                    residue = match_list[0][0] + match_list[0][2]
                    chain = match_list[0][1] + match_list[0][3]
                    if chain in chain_dict:
                        chain_dict[chain] += long2short[residue]
                    else:
                        chain_dict[chain] = long2short[residue]
                        chain_list.append(chain)

        with open(self.outfile, 'w') as fout:
            for chain in chain_list:
                fout.write(f'>{self.infile}:{chain}\n{chain_dict[chain]}\n')
