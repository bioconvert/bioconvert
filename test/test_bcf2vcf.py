from bioconvert.bcf2vcf import BCF2VCF
from easydev import TempFile, md5
import pytest

from . import test_dir


@pytest.mark.parametrize("method", BCF2VCF.available_methods)
def test_conv(method):
    infile = f"{test_dir}/data/bcf/test_bcf2vcf_v1.bcf"
    with TempFile(suffix=".vcf") as tempfile:
        convert = BCF2VCF(infile, tempfile.name)
        convert(method=method)
