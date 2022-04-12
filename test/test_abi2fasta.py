from bioconvert.abi2fasta import ABI2FASTA
import pytest
from easydev import TempFile, md5

from . import test_dir


@pytest.mark.parametrize("method", ABI2FASTA.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/abi/310.ab1"
    with TempFile(suffix=".fasta") as tempfile:
        convert = ABI2FASTA(infile, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) in ["547ab44aca8443b61e8ac2a5d0ce9a17"]
