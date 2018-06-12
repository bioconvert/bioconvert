import os
from bioconvert.bed2wiggle import BED2WIGGLE
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", BED2WIGGLE.available_methods)
def test_conv(method):
    infile = bioconvert_data("ucsc.bed")
    outfile = bioconvert_data("test_ucsc_bed2wiggle.wiggle")
    md5out = md5(outfile)

    with TempFile(suffix=".wiggle") as tempfile:
        convert = BED2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)

