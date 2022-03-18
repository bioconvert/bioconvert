import os
from bioconvert.bed2wiggle import BED2WIGGLE
from easydev import TempFile, md5
import pytest

from . import test_dir


@pytest.mark.parametrize("method", BED2WIGGLE.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/bed/ucsc.bed"
    outfile = f"{test_dir}/data/wiggle/test_ucsc_bed2wiggle.wiggle"
    md5out = md5(outfile)

    with TempFile(suffix=".wiggle") as tempfile:
        convert = BED2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)
