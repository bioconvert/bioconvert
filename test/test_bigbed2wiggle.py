import os
from bioconvert.bigbed2wiggle import BIGBED2WIGGLE
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", BIGBED2WIGGLE.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/bigbed/ucsc.bigbed"
    outfile = f"{test_dir}/data/wiggle/ucsc.wiggle"
    md5out = md5(outfile)


    with TempFile(suffix=".wiggle") as tempfile:
        convert = BIGBED2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)