import os
from bioconvert.bedgraph2wiggle import BEDGRAPH2WIGGLE
from bioconvert import TempFile, md5
import pytest

from . import test_dir


@pytest.mark.parametrize("method", BEDGRAPH2WIGGLE.available_methods)
def test_conv(method):
    infile = (
        f"{test_dir}/data/bedgraph/ucsc.bg"  # must have bg extension for wiggletools
    )
    outfile = f"{test_dir}/data/wiggle/test_bedgraph2wiggle.wiggle"
    md5out = md5(outfile)

    with TempFile(suffix=".wiggle") as tempfile:
        convert = BEDGRAPH2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)
