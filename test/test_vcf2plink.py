import os
import pytest
import tempfile
from bioconvert.vcf2plink import VCF2PLINK
from easydev import md5

from . import test_dir

@pytest.mark.parametrize("method", VCF2PLINK.available_methods)
def test_vcf2plink(method):
    infile = f"{test_dir}/data/vcf/plink_toy.vcf"
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink")
        converter = VCF2PLINK(infile, outfile)
        converter(method=method)

        assert os.path.isfile(outfile+".map")
        assert os.path.isfile(outfile+".ped")
        assert os.path.isfile(outfile+".log")