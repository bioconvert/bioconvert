import os
import pytest
import tempfile
from bioconvert.plink2vcf import PLINK2VCF
from easydev import md5

from . import test_dir


@pytest.mark.parametrize("method", PLINK2VCF.available_methods)
def test_plink2vcf(method):
    infile = os.path.splitext(f"{test_dir}/data/plink/plink_toy.ped")[0]
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink.vcf")
        converter = PLINK2VCF(infile, outfile)
        converter(method=method)

        assert os.path.isfile(outfile)
