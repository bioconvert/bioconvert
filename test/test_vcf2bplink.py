import os
import pytest
import tempfile
from bioconvert.vcf2bplink import VCF2BPLINK
from bioconvert import bioconvert_data
from easydev import md5


@pytest.mark.parametrize("method", VCF2BPLINK.available_methods)
def test_vcf2bplink(method):
    infile = bioconvert_data("plink_toy.vcf")
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink")
        converter = VCF2BPLINK(infile, outfile)
        converter(method=method)

        assert os.path.isfile(outfile+".bed")
        assert os.path.isfile(outfile+".bim")
        assert os.path.isfile(outfile+".fam")
