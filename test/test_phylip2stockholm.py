import pytest
from easydev import TempFile, md5

from bioconvert.phylip2stockholm import PHYLIP2STOCKHOLM

from . import test_dir


def test_phylip2stockholm_biopython():
    infile = f"{test_dir}/data/phylip/biopython.phylip"
    outfile = f"{test_dir}/data/stockholm/biopython.stockholm"
    with TempFile(suffix=".fasta") as tempfile:
        converter = PHYLIP2STOCKHOLM(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)


@pytest.mark.skipif(
    PHYLIP2STOCKHOLM._method_squizz.is_disabled, reason="missing dependencies"
)
def test_phylip2stockholm_squizz():
    infile = f"{test_dir}/data/phylip/squizz.phylip"
    outfile = f"{test_dir}/data/stockholm/squizz.stockholm"
    with TempFile(suffix=".stockholm") as tempfile:
        converter = PHYLIP2STOCKHOLM(infile, tempfile.name)
        converter(method="squizz")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(outfile)
