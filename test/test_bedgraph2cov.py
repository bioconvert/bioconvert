import pytest
from easydev import TempFile, md5

from bioconvert.bedgraph2cov import BEDGRAPH2COV

from . import test_dir

@pytest.mark.parametrize("method", BEDGRAPH2COV.available_methods)
def test_bedgraph2cov(method):
    infile = f"{test_dir}/data/bedgraph/test_bedgraph2bed.bedgraph"
    with TempFile(suffix=".cov") as tempfile:
        converter = BEDGRAPH2COV(infile, tempfile.name)
        converter(method=method)

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == "a8cc8b0fd2f2fd028424dc8969a0b8b6"