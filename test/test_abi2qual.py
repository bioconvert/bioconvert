from bioconvert.abi2qual import ABI2QUAL
from bioconvert import bioconvert_data
import pytest
from easydev import TempFile, md5


@pytest.mark.parametrize("method", ABI2QUAL.available_methods)
def test_conv(method):
    infile = bioconvert_data("ABI/310.ab1")
    with TempFile(suffix=".fastq") as tempfile:
        convert = ABI2QUAL(infile, tempfile.name)
        convert()
        assert md5(tempfile.name) in ["bdb0bb2dd00042733f145ba918dbab08"]
