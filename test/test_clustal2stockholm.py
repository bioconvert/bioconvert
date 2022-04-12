import os
import pytest
from easydev import TempFile, md5

from bioconvert.clustal2stockholm import CLUSTAL2STOCKHOLM

from . import test_dir


def test_clustal2stockholm_biopython():
    infile = f"{test_dir}/data/clustal/biopython.clustal"
    outfile = f"{test_dir}/data/stockholm/biopython.stockholm"
    with TempFile(suffix=".stockholm") as tempfile:
        converter = CLUSTAL2STOCKHOLM(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    CLUSTAL2STOCKHOLM._method_squizz.is_disabled, reason="missing dependencies"
)
def test_clustal2stockholm_squizz():
    infile = f"{test_dir}/data/clustal/squizz.clustal"
    outfile = f"{test_dir}/data/stockholm/squizz.stockholm"
    with TempFile(suffix=".stockholm") as tempfile:
        converter = CLUSTAL2STOCKHOLM(infile, tempfile.name)
        converter(method="squizz")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
