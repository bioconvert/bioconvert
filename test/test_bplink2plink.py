import os
import pytest
import tempfile
from bioconvert.bplink2plink import BPLINK2PLINK
from bioconvert import bioconvert_data
from easydev import md5


@pytest.mark.parametrize("method", BPLINK2PLINK.available_methods)
def test_bplink2plink(method):
    infile = os.path.splitext(bioconvert_data("plink_toy.bed"))[0]
    expected_outfile = os.path.splitext(bioconvert_data("plink_toy.ped"))[0]
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink")
        converter = BPLINK2PLINK(infile, outfile)
        converter(method=method)

        assert md5(outfile + ".ped") == md5(expected_outfile + ".ped")
        assert md5(outfile + ".map") == md5(expected_outfile + ".map")
