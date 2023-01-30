import pytest
from bioconvert import TempFile, md5

from bioconvert.twobit2fasta import TWOBIT2FASTA

from . import test_dir


@pytest.mark.parametrize("method", TWOBIT2FASTA.available_methods)
def test_twobit2fasta_ucsc(method):
    infile = f"{test_dir}/data/twobit/ucsc.2bit"
    outfile = f"{test_dir}/data/fasta/ucsc.fasta"
    with TempFile(suffix=".fasta") as tempfile:
        converter = TWOBIT2FASTA(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
