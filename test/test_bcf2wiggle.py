import os
from bioconvert.bcf2wiggle import BCF2WIGGLE
from bioconvert import TempFile, md5
import pytest

from . import test_dir


@pytest.mark.parametrize("method", BCF2WIGGLE.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/bcf/test_bcf2vcf_v1.bcf"
    outfile = f"{test_dir}/data/wiggle/test_bcf2wiggle.wiggle"
    md5out = md5(outfile)

    with TempFile(suffix=".wiggle") as tempfile:
        convert = BCF2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)
