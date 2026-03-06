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
"""Convert :term:`MOL2` to :term:`SMILES` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

__all__ = ["MOL22SMILES"]


class MOL22SMILES(ConvBase):
    """Convert :term:`MOL2` file to :term:`SMILES` file

    Convert Tripos MOL2 small molecule structure file into SMILES format.

    Methods available are based on rdkit [RDKIT]_ or openbabel [OPENBABEL]_.
    """

    #: Default value
    _default_method = "rdkit"

    def __init__(self, infile, outfile, *args, **kargs):
        """
        :param str infile: The path to the input MOL2 file.
        :param str outfile: The path to the output SMILES file.
        """
        super(MOL22SMILES, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="rdkit")
    def _method_rdkit(self, *args, **kwargs):
        """Convert MOL2 to SMILES using the rdkit library.

        `RDKit Documentation <https://www.rdkit.org/docs/>`_"""
        from rdkit import Chem

        mol = Chem.MolFromMol2File(self.infile, sanitize=True, removeHs=True)
        if mol is None:
            raise ValueError(f"Could not read molecule from {self.infile}")
        smiles = Chem.MolToSmiles(mol)
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
        with open(self.outfile, "w") as fout:
            fout.write(f"{smiles}\t{name}".rstrip() + "\n")

    @requires(external_binary="obabel")
    def _method_openbabel(self, *args, **kwargs):
        """Convert MOL2 to SMILES using the openbabel command-line tool.

        `OpenBabel Documentation <http://openbabel.org/>`_"""
        cmd = f"obabel {self.infile} -osmi -O {self.outfile}"
        self.execute(cmd)
