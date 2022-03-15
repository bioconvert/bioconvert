from bioconvert.abi2qual import ABI2QUAL
import pytest
from easydev import TempFile, md5

from . import test_dir

@pytest.mark.parametrize("method", ABI2QUAL.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/abi/310.ab1"
    with TempFile(suffix=".fastq") as tempfile:
        convert = ABI2QUAL(infile, tempfile.name)
        convert()
        assert md5(tempfile.name) in ["bdb0bb2dd00042733f145ba918dbab08"]
