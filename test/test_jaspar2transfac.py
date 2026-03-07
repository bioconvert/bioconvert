import pytest

from bioconvert.jaspar2transfac import JASPAR2TRANSFAC
from bioconvert import TempFile, md5

from . import test_dir


@pytest.mark.skipif(
    JASPAR2TRANSFAC._method_biopython.is_disabled, reason="missing dependencies"
)
def test_jaspar2transfac_biopython():
    infile = f"{test_dir}/data/jaspar/test_motifs.jaspar"
    expected_outfile = f"{test_dir}/data/transfac/test_motifs_from_jaspar.transfac"
    with TempFile(suffix=".transfac") as tempfile:
        converter = JASPAR2TRANSFAC(infile, tempfile.name)
        converter(method="biopython")

        # Check that the output is correct with a checksum
        assert md5(tempfile.name) == md5(expected_outfile)
