import os
import pytest
import tempfile
from bioconvert.bplink2vcf import BPLINK2VCF
from bioconvert import bioconvert_data
from easydev import md5


@pytest.mark.parametrize("method", BPLINK2VCF.available_methods)
def test_bplink2vcf(method):
    infile = os.path.splitext(bioconvert_data("plink_toy.bed"))[0]
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink.vcf")
        converter = BPLINK2VCF(infile, outfile)
        converter(method=method)

        assert os.path.isfile(outfile)
