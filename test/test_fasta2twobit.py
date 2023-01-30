import pytest
from bioconvert import TempFile, md5

from bioconvert.fasta2twobit import FASTA2TWOBIT

from . import test_dir


@pytest.mark.parametrize("method", FASTA2TWOBIT.available_methods)
def test_fasta2twobit_ucsc(method):
    infile = f"{test_dir}/data/fasta/ucsc.fasta"
    outfile = f"{test_dir}/data/twobit/ucsc.2bit"
    with TempFile(suffix=".2bit") as tempfile:
        converter = FASTA2TWOBIT(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
