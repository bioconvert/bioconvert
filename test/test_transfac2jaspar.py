import pytest

from bioconvert.transfac2jaspar import TRANSFAC2JASPAR
from bioconvert import TempFile, md5

from . import test_dir


@pytest.mark.skipif(
    TRANSFAC2JASPAR._method_biopython.is_disabled, reason="missing dependencies"
)
def test_transfac2jaspar_biopython():
    infile = f"{test_dir}/data/transfac/test_motifs.transfac"
    expected_outfile = f"{test_dir}/data/jaspar/test_motifs.jaspar"
    with TempFile(suffix=".jaspar") as tempfile:
        converter = TRANSFAC2JASPAR(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outfile)
