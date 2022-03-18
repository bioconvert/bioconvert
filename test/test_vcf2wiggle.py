from bioconvert.vcf2wiggle import VCF2WIGGLE
from easydev import TempFile, md5
import pytest

from . import test_dir

@pytest.mark.parametrize("method", VCF2WIGGLE.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/vcf/test_vcf2bcf_v1.vcf"
    outfile = f"{test_dir}/data/wiggle/test_vcf2wiggle.wiggle"
    md5out = md5(outfile)


    with TempFile(suffix=".wiggle") as tempfile:
        convert = VCF2WIGGLE(infile, tempfile.name)
        convert(method=method)

        assert md5(tempfile.name) == md5out, "{} failed".format(method)