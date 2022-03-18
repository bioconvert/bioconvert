import pytest

from bioconvert.vcf2bcf import VCF2BCF
from easydev import TempFile, md5

from . import test_dir


@pytest.mark.skipif(VCF2BCF._method_bcftools.is_disabled, reason="missing dependencies")
def test_conv():
    infile = f"{test_dir}/data/vcf/test_vcf2bcf_v1.vcf"
    with TempFile(suffix=".bcf") as tempfile:
        convert = VCF2BCF(infile, tempfile.name)
        convert(method="bcftools")
