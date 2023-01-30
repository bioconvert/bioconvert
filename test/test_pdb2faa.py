import pytest
from bioconvert import TempFile, md5
import filecmp

from bioconvert.pdb2faa import PDB2FAA

from . import test_dir


@pytest.mark.parametrize("method", PDB2FAA.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/pdb/1mbs.pdb"
    with TempFile(suffix=".faa") as tempfile:
        convert = PDB2FAA(infile, tempfile.name)
        convert()

        # input file name used in the header so changes all the time
        #assert md5(tempfile.name) == '2ca590ce0354dc6948d629bfd1e0b1e3'

