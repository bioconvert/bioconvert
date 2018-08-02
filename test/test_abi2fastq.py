from bioconvert.abi2fastq import ABI2FASTQ
from bioconvert import bioconvert_data
import pytest
from easydev import TempFile, md5


@pytest.mark.parametrize("method", ABI2FASTQ.available_methods)
def test_conv(method):
    infile = bioconvert_data("ABI/310.ab1")
    with TempFile(suffix=".fastq") as tempfile:
        convert = ABI2FASTQ(infile, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        # Note that awk wrap fastq while python method does not 
        assert md5(tempfile.name) in \
            [ "a1c5028da7c0429fa5d9e8b6ef9d3691"]
 
