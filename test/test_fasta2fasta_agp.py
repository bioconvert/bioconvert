from bioconvert.fasta2fasta_agp import FASTA2FASTA_AGP
from bioconvert import TempFile, md5
import pytest

from . import test_dir

@pytest.mark.parametrize("method", FASTA2FASTA_AGP.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/fasta/biopython.fasta"

    with TempFile(suffix=".fasta") as out1, TempFile(suffix=".agp") as out2:
        convert = FASTA2FASTA_AGP(infile, (out1.name, out2.name))
        convert(method=method)

