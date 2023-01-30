import os
import pytest
import tempfile
from bioconvert.bplink2plink import BPLINK2PLINK
from bioconvert import md5

from . import test_dir


@pytest.mark.parametrize("method", BPLINK2PLINK.available_methods)
def test_bplink2plink(method):
    infile = os.path.splitext(f"{test_dir}/data/bed/plink_toy.bed")[0]
    expected_outfile = os.path.splitext(f"{test_dir}/data/plink/plink_toy.ped")[0]
    with tempfile.TemporaryDirectory() as tmpdirname:
        outfile = os.path.join(tmpdirname, "plink")
        converter = BPLINK2PLINK(infile, outfile)
        converter(method=method)

        assert md5(outfile + ".ped") == md5(expected_outfile + ".ped")
        assert md5(outfile + ".map") == md5(expected_outfile + ".map")
