import os
import pytest
from easydev import TempFile, md5

from bioconvert.stockholm2clustal import STOCKHOLM2CLUSTAL

from . import test_dir


def test_stockholm2clustal_biopython():
    infile = f"{test_dir}/data/stockholm/biopython.stockholm"
    outfile = f"{test_dir}/data/clustal/biopython.clustal"
    with TempFile(suffix=".clustal") as tempfile:
        converter = STOCKHOLM2CLUSTAL(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    STOCKHOLM2CLUSTAL._method_squizz.is_disabled, reason="missing dependencies"
)
def test_stockholm2clustal_squizz():
    infile = f"{test_dir}/data/stockholm/squizz.stockholm"
    outfile = f"{test_dir}/data/clustal/squizz.clustal"
    with TempFile(suffix=".clustal") as tempfile:
        converter = STOCKHOLM2CLUSTAL(infile, tempfile.name)
        converter(method="squizz")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
