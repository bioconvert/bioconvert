import os
import pytest
import tempfile
from bioconvert.plink2bplink import PLINK2BPLINK
from bioconvert import bioconvert_data
from easydev import md5


@pytest.mark.parametrize("method", PLINK2BPLINK.available_methods)
def test_plink2bplink(method):
    infile = os.path.splitext(bioconvert_data("plink_toy.ped"))[0]
    expected_outfile = os.path.splitext(bioconvert_data("plink_toy.bed"))[0]
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink")
        converter = PLINK2BPLINK(infile, outfile)
        converter(method=method)

        assert md5(outfile + ".bed") == md5(expected_outfile + ".bed")
        assert md5(outfile + ".bim") == md5(expected_outfile + ".bim")
        assert md5(outfile + ".fam") == md5(expected_outfile + ".fam")
