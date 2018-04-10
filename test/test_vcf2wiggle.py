import os
from bioconvert.vcf2wiggle import VCF2WIGGLE
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest



@pytest.mark.parametrize("method", VCF2WIGGLE.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_vcf2bcf_v1.vcf")
    outfile = bioconvert_data("test_vcf2bcf_v1.wiggle")
    md5out = md5(outfile)


    with TempFile(suffix=".wiggle") as tempfile:
        convert = VCF2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)

