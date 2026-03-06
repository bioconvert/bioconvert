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
"""Convert :term:`SDF` to :term:`SMILES` format"""
from bioconvert import ConvBase
from bioconvert.core.decorators import requires

__all__ = ["SDF2SMILES"]


class SDF2SMILES(ConvBase):
    """Convert :term:`SDF` file to :term:`SMILES` file

    Convert MDL SDF (Structure-Data File) small molecule structure file
    into SMILES format. Each molecule in the SDF file is written as one
    SMILES entry per line.

    Methods available are based on rdkit [RDKIT]_ or openbabel [OPENBABEL]_.
    """

    #: Default value
    _default_method = "rdkit"

    def __init__(self, infile, outfile, *args, **kargs):
        """
        :param str infile: The path to the input SDF file.
        :param str outfile: The path to the output SMILES file.
        """
        super(SDF2SMILES, self).__init__(infile, outfile, *args, **kargs)

    @requires(python_library="rdkit")
    def _method_rdkit(self, *args, **kwargs):
        """Convert SDF to SMILES using the rdkit library.

        `RDKit Documentation <https://www.rdkit.org/docs/>`_"""
        from rdkit import Chem

        suppl = Chem.SDMolSupplier(self.infile, sanitize=True, removeHs=True)
        with open(self.outfile, "w") as fout:
            for mol in suppl:
                if mol is None:
                    continue
                smiles = Chem.MolToSmiles(mol)
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else ""
                line = f"{smiles}\t{name}".rstrip()
                fout.write(line + "\n")

    @requires(external_binary="obabel")
    def _method_openbabel(self, *args, **kwargs):
        """Convert SDF to SMILES using the openbabel command-line tool.

        `OpenBabel Documentation <http://openbabel.org/>`_"""
        cmd = f"obabel {self.infile} -osmi -O {self.outfile}"
        self.execute(cmd)
