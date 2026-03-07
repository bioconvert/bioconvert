import pytest
from bioconvert import TempFile

from bioconvert.sdf2smiles import SDF2SMILES

from . import test_dir


@pytest.mark.parametrize("method", SDF2SMILES.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/sdf/aspirin.sdf"
    with TempFile(suffix=".smi") as tempfile:
        convert = SDF2SMILES(infile, tempfile.name)
        convert(method=method)

        with open(tempfile.name) as f:
            content = f.read().strip()

        assert content, "Output SMILES file should not be empty"
        # The first field (before tab or end of line) should be a valid SMILES
        smiles = content.split("\t")[0].split()[0]
        assert smiles, "SMILES string should not be empty"
        # Check for expected atoms in aspirin
        assert "O" in smiles, "Aspirin SMILES should contain oxygen"
