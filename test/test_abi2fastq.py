from bioconvert.abi2fastq import ABI2FASTQ
import pytest
from easydev import TempFile, md5

from . import test_dir


@pytest.mark.parametrize("method", ABI2FASTQ.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/abi/310.ab1"

    with TempFile(suffix=".fastq") as tempfile:
        convert = ABI2FASTQ(infile, tempfile.name)
        convert()
        # Check that the output is correct with a checksum
        assert md5(tempfile.name) in ["a1c5028da7c0429fa5d9e8b6ef9d3691"]
