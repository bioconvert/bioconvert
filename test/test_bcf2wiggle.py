import os
from bioconvert.bam2wiggle import BCF2WIGGLE
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", BCF2WIGGLE.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_bcf2vcf_v1.bcf")
    outfile = bioconvert_data("test_bam2wiggle.wiggle")
    md5out = md5(outfile)


    with TempFile(suffix=".wiggle") as tempfile:
        convert = BCF2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)

