import os
from bioconvert.bam2wiggle import BAM2WIGGLE
from easydev import TempFile, md5
import pytest

from . import test_dir

@pytest.mark.parametrize("method", BAM2WIGGLE.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/bam/test_measles.sorted.bam"
    outfile = f"{test_dir}/data/wiggle/test_bam2wiggle.wiggle"
    md5out = md5(outfile)

    with TempFile(suffix=".wiggle") as tempfile:
        convert = BAM2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)

