import os
import pytest
import tempfile
from bioconvert.bplink2vcf import BPLINK2VCF
from easydev import md5

from . import test_dir

@pytest.mark.parametrize("method", BPLINK2VCF.available_methods)
def test_bplink2vcf(method):
    infile = os.path.splitext(f"{test_dir}/data/bed/plink_toy.bed")[0]
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink.vcf")
        converter = BPLINK2VCF(infile, outfile)
        converter(method=method)

        assert os.path.isfile(outfile)