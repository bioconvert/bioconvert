from bioconvert.abi2fasta import ABI2FASTA
from bioconvert import bioconvert_data
import pytest
from easydev import TempFile, md5


@pytest.mark.parametrize("method", ABI2FASTA.available_methods)
def test_conv(method):
    infile = bioconvert_data("ABI/310.ab1")
    with TempFile(suffix=".fasta") as tempfile:
        convert = ABI2FASTA(infile, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) in \
            ["547ab44aca8443b61e8ac2a5d0ce9a17"]
