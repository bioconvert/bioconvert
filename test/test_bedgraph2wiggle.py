import os
from bioconvert.bam2wiggle import BEDGRAPH2WIGGLE
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", BEDGRAPH2WIGGLE.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_measles.bedgraph")
    outfile = bioconvert_data("test_bedgraph2wiggle.wiggle")
    md5out = md5(outfile)


    with TempFile(suffix=".wiggle") as tempfile:
        convert = BEDGRAPH2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)

