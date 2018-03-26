import pytest

from bioconvert.vcf2bcf import VCF2BCF
from bioconvert import bioconvert_data
from easydev import TempFile, md5


@pytest.mark.skipif(VCF2BCF._method_bcftools.is_disabled, reason="missing dependencies")
def test_conv():
    infile = bioconvert_data("test_vcf2bcf_v1.vcf")
    with TempFile(suffix=".bcf") as tempfile:
        convert = VCF2BCF(infile, tempfile.name)
        convert(method="bcftools")
