from bioconvert.bcf2vcf import BCF2VCF
from bioconvert import bioconvert_data
from easydev import TempFile, md5
import pytest

@pytest.mark.parametrize("method", BCF2VCF.available_methods)
def test_conv(method):
    infile = bioconvert_data("test_bcf2vcf_v1.bcf")
    with TempFile(suffix=".vcf") as tempfile:
        convert = BCF2VCF(infile, tempfile.name)
        convert(method=method)


